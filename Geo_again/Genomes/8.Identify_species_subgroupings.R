#!/usr/local/bin/Rscript
library(parallel)
library(stringr)
library(dplyr)

bullShitTableRead <- function(fileName, header = NA) {
	# Read in file
	con 	<- file(fileName)
	data	<- readLines(con = con)
	close(con)

	# Sanitise eol
	eolFixed 	<- str_replace(data, "(\t{1}[|]{1}$)", "")

	# Split by stupid separator + combine to df
	records		<- str_split(eolFixed, "(\t{1}[|]{1}\t{1})", simplify = TRUE)
	records_df	<- as.data.frame(records, stringsAsFactors = FALSE)

	# Name if want
	if(!is.na(header) && ncol(records_df) == length(header)) {
		names(records_df) <- header
	}

	return(records_df)
}

# ----------------------------------------------------------------------------------- #

direct		<- "/Users/aesin/Desktop/Geo_again"
genome_dir	<- file.path(direct, "Genomes")
taxdmp_dir	<- file.path(genome_dir, "taxdmp")

## Taxdmp files
names_file	<- file.path(taxdmp_dir, "names.dmp")
nodes_file 	<- file.path(taxdmp_dir, "nodes.dmp")

## Translation tables
accass_tax_file	<- file.path(genome_dir, "Genome_lists", "Acc_ass_taxid_table.tsv")
binom_tax_file	<- file.path(genome_dir, "Genome_lists", "Taxid_refinedBinomial_table.tsv")

## Header names for taxdmp file
namesHeader	<- c("taxid", "name", "uniqueName", "nameClass")
nodesHeader	<- c("taxid", "parentTaxid", "nodeRank", "EMBLcode", "divisionID", "divFlag", "genCodeID", "genCodeFlag", "mitoGenCode", "mitoFlag", "genbankHidFlag", "hiddenSubTreeFlag", "comments")

## Read in data
names_data		<- bullShitTableRead(fileName = names_file, header = namesHeader)
nodes_data		<- bullShitTableRead(fileName = nodes_file, header = nodesHeader)
accass_tax_tbl	<- read.table(accass_tax_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
names(accass_tax_tbl)	<- c("AccAss", "Taxid")
binom_tax_tbl	<- read.table(binom_tax_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# ----------------------------------------------------------------------------------- #

taxid_list	<- accass_tax_tbl$Taxid
taxidNum	<- length(taxid_list)

numCores	<- 20
clust		<- makeCluster(numCores, type = "FORK")
perTaxidSpecies	<- parLapply(clust, 1:length(taxid_list), function(taxidIndex) {

	message(paste0("\rProcessing taxid relatedness... ", taxidIndex, " // ", taxidNum), appendLF = FALSE)
	
	taxid		<- taxid_list[taxidIndex]
	childTaxid	<- taxid


	while(TRUE) {
		parentTaxid	<- nodes_data$parentTaxid[which(nodes_data$taxid == childTaxid)]
		childTaxid	<- ifelse(length(which(taxid_list == parentTaxid)) > 0, parentTaxid, childTaxid)

		if (!identical(parentTaxid, childTaxid)) {
			break
		}
	}

	# 3 taxids are missing from the binomial table - account for this
	binomIndex	<- which(binom_tax_tbl$Taxid == childTaxid)
	parentBinom	<- ifelse(length(binomIndex) == 1, binom_tax_tbl$Binomial[binomIndex], NA)

	# Output row
	out_df		<- data.frame(taxid = taxid, internalParent = as.numeric(childTaxid), parentName = parentBinom, stringsAsFactors = FALSE)
	return(out_df)
})

perTaxidSpecies_df	<- bind_rows(perTaxidSpecies)
stopCluster(clust)

write.table(x = perTaxidSpecies_df, file = file.path(genome_dir, "Genome_lists", "Species_groupings.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

