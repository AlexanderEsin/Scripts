#!/usr/bin/env Rscript
library(dplyr)
library(RSQLite)

# ------------------------------------------------------------------------------------- #
genomeRelativePosition_format	<- function(row, genePosName) {

	oriStrand	<- row[["dnaA_strand"]]
	geneStart	<- as.numeric(row[[genePosName]])
	genome_l	<- as.numeric(row[["GenomeLen"]])

	if (is.na(geneStart)) {
		return(NA)
	}

	if (identical(oriStrand, "Forward")) {

		oriStart	<- as.numeric(row[["dnaA_start"]])

		# Get the position relative to the origin
		relGenePosition	<- geneStart - oriStart

		# If the relative location is negative then it lies upstream of the origin
		if (relGenePosition < 0) {
			relGenePosition	<- relGenePosition + genome_l
		}
	} else if (identical(oriStrand, "Reverse")) {

		oriEnd	<- as.numeric(row[["dnaA_end"]])

		# Get the position relative to the origin
		relGenePosition	<- geneStart - oriEnd

		# If it's negative it lies downstream of the origin. Take the absolute value and leave as is #
		# If it's positive it lies upstream. We subtract it's position from the genome size. This will result in chromStart being > chromEnd #
		relGenePosition	<- ifelse(relGenePosition < 0, abs(relGenePosition), genome_l - relGenePosition)
	} else {
		stop("Origin strand information is missing")
	}

	## Get the fractional position
	fractGenePosition	<- round(relGenePosition / genome_l, digits = 8)
	return(fractGenePosition)
}

# ------------------------------------------------------------------------------------- #

master_path			<- "/Users/aesin/Desktop/Geo_again"
genome_path			<- file.path(master_path, "Genomes")
genomeList_path		<- file.path(genome_path, "Genome_lists")
islandView_path		<- file.path(genome_path, "islandViewer")

giProcess_path		<- file.path(master_path, "HGT_position", "GI_data")
if (!dir.exists(giProcess_path)) dir.create(giProcess_path)

allProtDB_path		<- file.path(master_path, "All_prot_db_new")
# ------------------------------------------------------------------------------------- #

# Get the AccAss-to-Chromosome table
accAssToChrom_file	<- file.path(genomeList_path, "AG_acc_ass_chromosome.txt")
accAssToChrom_data	<- read.table(accAssToChrom_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Get the AccAss-to-Taxid table
accAssToTaxid_file	<- file.path(genomeList_path, "Acc_ass_taxid_table.tsv")
accAssToTaxid_data	<- read.table(accAssToTaxid_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
names(accAssToTaxid_data)	<- c("AccAss", "Taxid")

# Get the Taxid-to-dnaA table
taxidToDnaA_file	<- file.path(genomeList_path, "AG_dnaA_startEnd_positions.tsv")
taxidToDnaA_data	<- read.table(taxidToDnaA_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
names(taxidToDnaA_data)[2:4]	<- c("dnaA_start", "dnaA_end", "dnaA_strand")

# Combine to have Chromosome linked to Taxid
combined_df			<- left_join(accAssToChrom_data, accAssToTaxid_data, by = "AccAss")
combined_df			<- left_join(combined_df, taxidToDnaA_data, by = "Taxid")

# ------------------------------------------------------------------------------------- #

# Open All_prot database
conn			<- dbConnect(RSQLite::SQLite(), allProtDB_path)

# Extract taxid and genome_l
genomeLen_tbl		<- dbSendQuery(conn, 'SELECT taxid, binomial, genome_l FROM t1 WHERE taxid = :taxids LIMIT 1')
dbBind(genomeLen_tbl, param = list(taxids = combined_df$Taxid))
genomeLen_df		<- dbFetch(genomeLen_tbl)
dbClearResult(genomeLen_tbl)
names(genomeLen_df)	<- c("Taxid", "Binomial", "GenomeLen")

# Add length to combined genome
combined_df			<- left_join(combined_df, genomeLen_df, by = "Taxid")

# Disconnect
dbDisconnect(conn)

# ------------------------------------------------------------------------------------- #

# Open the islandViewer data
islandView_file		<- file.path(islandView_path, "all_gis_islandviewer_iv4.txt")
islandView_data		<- read.table(islandView_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
names(islandView_data)[1]	<- "Chromosome"

# Combine island view data to AG accAss etc...
allCombined_df		<- left_join(combined_df, islandView_data, by = "Chromosome")

# Get the GI relative start and end positions
GI_relStart			<- apply(allCombined_df, 1, genomeRelativePosition_format, genePosName = "Start")
GI_relEnd			<- apply(allCombined_df, 1, genomeRelativePosition_format, genePosName = "End")

# Because in some cases, the Ori strand is reverse, GI start and end can be backwards.
# We'll just take the lower relative position as the start
adjGI_boundaries	<- lapply(1:nrow(allCombined_df), function(entryIndex) {
	relStart	<- GI_relStart[entryIndex]
	relEnd		<- GI_relEnd[entryIndex]

	if (is.na(relStart)) {
		out_df	<- data.frame(GI_relStart = NA, GI_relEnd = NA, stringsAsFactors = FALSE)
		return(out_df)
	}

	adjStart	<- ifelse(relStart < relEnd, relStart, relEnd)
	adjEnd		<- ifelse(relEnd > relStart, relEnd, relStart)

	out_df	<- data.frame(GI_relStart = adjStart, GI_relEnd = adjEnd, stringsAsFactors = FALSE)
	return(out_df)
})
adjGI_boundaries_df	<- bind_rows(adjGI_boundaries)

# Bind rows
allCombined_withRelGI_df	<- cbind(allCombined_df, adjGI_boundaries_df)

# Write out the table
write.table(allCombined_withRelGI_df, file = file.path(giProcess_path, "perGenome_relGIBoundaries.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)



# ------------------------------------------------------------------------------------- #









