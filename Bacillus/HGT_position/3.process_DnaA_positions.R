#!/usr/bin/env Rscript

# Load master variables and HGT position functions
invisible(sapply(HGTPos.bac, source, .GlobalEnv))

require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load(tidyverse, fs, RSQLite)

master_path			<- "/Users/aesin/Desktop/Bacillus/"
genome_list_path	<- file.path(master_path, "Core_genomes", "Genome_lists")

dbConn			<- dbConnect(RSQLite::SQLite(), allProtDB_path)

# ------------------------------------ #
dnaA_locus_file	<- file.path(genome_list_path, "Core_origin_gene_list.tsv")
dnaA_locus_data	<- read_tsv(dnaA_locus_file)

core_dnaA_query	<- dbSendQuery(dbConn, 'SELECT protID, locus, OrthGroup FROM t1 WHERE locus == :locusList')
dbBind(core_dnaA_query, param = list(locusList = dnaA_locus_data$dnaATag))
core_dnaA_tbl	<- as_tibble(dbFetch(core_dnaA_query))
dbClearResult(core_dnaA_query)

# Extract the family / group to which DnaA proteins belong
dnaA_family		<- core_dnaA_tbl %>% distinct(OrthGroup) %>% pull(OrthGroup)

# ------------------------------------ #
# Get the protIDs from the base family (1e-10)
baseGroup_file	<- file.path(master_path, "Family_groups", "-10", "Final_groups.tsv")
baseGroup_data	<- read_tsv(baseGroup_file, col_names = FALSE)
all_dnaA_protID	<- baseGroup_data %>% filter(X1 == dnaA_family) %>% select(-1) %>% unlist(, use.names = FALSE)


# ------------------------------------ #
# Get dnaA data for all proteins in this family
all_dnaA_query	<- dbSendQuery(dbConn, 'SELECT * FROM t1 WHERE protID == :protID_list;')
dbBind(all_dnaA_query, param = list(protID_list = all_dnaA_protID))
all_dnaA_tbl	<- as_tibble(dbFetch(all_dnaA_query))
dbClearResult(all_dnaA_query)

# Clean DNAa table to contain only taxid entries with one DNAa gene which is not on a plasmid
dnaA_clean_tbl	<- all_dnaA_tbl %>%
	count(taxid) %>%
	filter(n == 1) %>%
	left_join(all_dnaA_tbl, by = "taxid") %>%
	filter(plasmid == "F")

dnaA_clean_trm	<- dnaA_clean_tbl %>%
	select(taxid, gene_start, gene_end, strand) %>%
	rename(oriStart = gene_start, oriEnd = gene_end, oriStrand = strand)

# Read in the species subgroupings file (result of Scripts/Geo_again/Genomes/8.Identify_species_subgroupings.R)
# speciesSubgroups_file	<- file.path(genome_path, "Genome_lists", "Species_groupings.tsv")
# speciesSubgroups_df		<- read.table(speciesSubgroups_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "")

# Save output
saveRDS(object = dnaA_clean_trm, file = file.path(positionData_path, "bySpecies_dnaA_data.rds"))
# Close db
dbDisconnect(dbConn)