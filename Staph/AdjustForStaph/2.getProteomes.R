#!/usr/bin/env Rscript

# Any libraries explicitly used in the script
require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load(tidyverse, RSQLite, fs)

# ----------------------------------------------------- #
master			<- file.path("/Users/aesin/Desktop/Staph")

# List of all proteomes
proteome_dir	<- file.path("/Users/aesin/Desktop/Geo_again/Proteomes/")
all_clean		<- file.path(proteome_dir, "Proteome_clean_fastas")
all_dup			<- file.path(proteome_dir, "Proteome_dup_fastas")

# List of genomes to remove
grp_list_dir	<- file.path(master, "Staph_genomes", "Genome_lists")
toRemove		<- read_tsv(file.path(grp_list_dir, "grpToRemove.tsv"))
accAss_tax_tbl	<- read_tsv(file.path(grp_list_dir, "grp_AccAssTaxid_table.tsv"))

# Create new directories for Clean and Duplicated proteomes
grpProteome_dir	<- file.path(master, "Proteomes")
clean_dir		<- file.path(grpProteome_dir, "Clean")
dup_dir			<- file.path(grpProteome_dir, "Dup")
dir_create(c(clean_dir, dup_dir))

# ----------------------------------------------------- #
# Copy all relevant clean and duplicated proteomes into the new directory
invisible(lapply(1:nrow(accAss_tax_tbl), function(row_index) {
	acc_ass		<- accAss_tax_tbl[row_index,] %>% pull(Acc_ass)
	file_name	<- paste0(acc_ass, "_proteome.fasta.gz")
	clean_file	<- file.path(all_clean, file_name)
	dup_file	<- file.path(all_dup, file_name)

	file_copy(c(clean_file, dup_file), c(clean_dir, dup_dir))

	return(NULL)
}))