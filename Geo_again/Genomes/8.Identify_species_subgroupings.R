#!/usr/local/bin/Rscript
library(parallel)
library(stringr)
library(dplyr)
library(data.table)

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
# Paths & files
direct		<- "/Users/aesin/Desktop/Geo_again"
genome_dir	<- file.path(direct, "Genomes")
taxdmp_dir	<- file.path(genome_dir, "taxdmp")

# Taxdmp files. If the RDS objects exist, use those:
if (file.exists(file.path(taxdmp_dir, "names.rds"))) {
	# Read in RDS objects
	names_dt	<- readRDS(file.path(taxdmp_dir, "names.rds"))
	nodes_dt	<- readRDS(file.path(taxdmp_dir, "nodes.rds"))
} else {
	message("\nNames and nodes RDS files not found: reading in from raw...", appendLF = FALSE)
	# Define file names
	names_file	<- file.path(taxdmp_dir, "names.dmp")
	nodes_file 	<- file.path(taxdmp_dir, "nodes.dmp")

	# Header names for taxdmp files
	namesHeader	<- c("taxid", "name", "uniqueName", "nameClass")
	nodesHeader	<- c("taxid", "parentTaxid", "nodeRank", "EMBLcode", "divisionID", "divFlag", "genCodeID", "genCodeFlag", "mitoGenCode", "mitoFlag", "genbankHidFlag", "hiddenSubTreeFlag", "comments")

	# Read in the files (takes a while)
	names_data		<- bullShitTableRead(fileName = names_file, header = namesHeader)
	nodes_data		<- bullShitTableRead(fileName = nodes_file, header = nodesHeader)

	# Change the relevant columns to numeric
	mode(names_data$taxid)			<- "numeric"
	mode(nodes_data$taxid) 			<- "numeric"
	mode(nodes_data$parentTaxid)	<- "numeric"

	# Convert to datatable
	names_dt	<- as.data.table(names_data)
	nodes_dt	<- as.data.table(nodes_data)

	# Set keys
	setkey(nodes_dt, taxid, parentTaxid)
	setkey(names_dt, taxid)

	# Write the names and nodes data tables to R objects to avoid reading in everytime
	saveRDS(object = names_dt, file = file.path(taxdmp_dir, "names.rds"))
	saveRDS(object = nodes_dt, file = file.path(taxdmp_dir, "nodes.rds"))

	# Clean
	rm(list = c("names_data", "nodes_data")); gc()

	message("\rNames and nodes RDS files not found: reading in from raw... done & saved as RDS")
}


# Translation table files
accass_tax_file	<- file.path(genome_dir, "Genome_lists", "Acc_ass_taxid_table.tsv")
binom_tax_file	<- file.path(genome_dir, "Genome_lists", "Taxid_refinedBinomial_table.tsv")

# Read in translation tables
accass_tax_tbl	<- read.table(accass_tax_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
names(accass_tax_tbl)	<- c("AccAss", "Taxid")
binom_tax_tbl	<- read.table(binom_tax_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)


# ----------------------------------------------------------------------------------- #

# List of taxids for which to find top-most species identifications
taxid_list	<- accass_tax_tbl$Taxid
taxidNum	<- length(taxid_list)

# For each taxid in our dataset, we identify the topmost (species) level taxid.
perTaxidSpecies <- mclapply(1:length(taxid_list), function(taxidIndex) {

	taxid		<- taxid_list[taxidIndex]
	childTaxid	<- taxid
	childRank	<- nodes_dt[.(childTaxid), nodeRank]

	if (identical(childRank, "species")) {
		return(data.frame(taxid = taxid, parentSpecies = childTaxid, stringsAsFactors = FALSE))
	}

	while(TRUE) {
		parentTaxidVal	<- nodes_dt[.(childTaxid), parentTaxid]
		parentRank		<- nodes_dt[.(parentTaxidVal), nodeRank]

		if (identical(parentRank, "species")) {
			break
		} else {
			childTaxid	<- parentTaxidVal
		}
	}

	return(data.frame(taxid = taxid, parentSpecies = parentTaxidVal, stringsAsFactors = FALSE))
}, mc.cores = 20)

# Combine to df and close cluster
perTaxidSpecies_df	<- bind_rows(perTaxidSpecies)

# How many unique species do we have?
uniqueParentTaxids	<- unique(perTaxidSpecies_df$parentSpecies)

# Keep only the scientific name entries in the names table
names_dataReduct	<- names_dt[.(nameClass) == "scientific name"),]

# Add the species names to the groups
perTaxidSpeciesNames_df	<- left_join(perTaxidSpecies_df, names_dataReduct[,1:2], by = c("parentSpecies" = "taxid"))
names(perTaxidSpeciesNames_df)[3]	<- "parentName"

# Write out results
write.table(x = perTaxidSpeciesNames_df, file = file.path(genome_dir, "Genome_lists", "Species_groupings.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

