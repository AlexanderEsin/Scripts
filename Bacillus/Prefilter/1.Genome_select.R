#!/usr/bin/env Rscript

# Load master variables and HGT position functions
invisible(sapply(HGTPos.all, source, .GlobalEnv))

# Any libraries explicitly used in the script
require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load(tidyverse, ape, phangorn)

## ----------------------------------------------- ##

master		<- "/Users/aesin/Desktop/Bacillus"
bacillac	<- file.path(master, "Prefilter", "Bacillaceae")

# Translate binomials to taxids
taxid_binomd_tbl	<- read_tsv(file.path(master, "Genomes", "Genome_lists", "Taxid_refinedBinomial_table.tsv"))
acc_ass_taxid_tbl	<- read_tsv(file.path(master, "Genomes", "Genome_lists", "Acc_ass_taxid_table.tsv"), col_names = c("Acc_ass", "Taxid"))

tree_file	<- file.path(bacillac, 	"RAxML_binomial_SupTree_-50.txt")
tree_data	<- read.tree(tree_file)

# Selet the group comprising target group + outgroup
edge_tips	<- c("Bacillus_subtilis", "Bacillus_cohnii")
mrca_node	<- getMRCA(tree_data, edge_tips)

# Extract the clade we want and get the taxids
get_clade	<- extract.clade(tree_data, node = mrca_node)
get_taxids	<- taxid_binomd_tbl %>% filter(Binomial %in% get_clade$tip.label)
get_accAss	<- acc_ass_taxid_tbl %>% filter(Taxid %in% get_taxids$Taxid)

## ----------------------------------------------- ##

# Create output master dir
output_dir	<- file.path(master, "Prefilter", "Bacillus_withOut", "Proteomes")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

write_delim(get_accAss %>% select(Acc_ass), path = file.path(output_dir, "..", "Acc_ass_list.txt"))

# Transfer the relevant proteomes (clean and dup)
clean_dir	<- file.path(output_dir, "Clean")
dup_dir		<- file.path(output_dir, "Dup")
invisible(lapply(c(clean_dir, dup_dir), function(x) if(!dir.exists(x)) dir.create(x)))

invisible(lapply(get_accAss$Acc_ass, function(acc_ass) {

	prot_name	<- paste0(acc_ass, "_proteome.fasta.gz")
	clean_prot	<- file.path(bacillac, "Proteomes", "Clean", prot_name)
	dup_prot	<- file.path(bacillac, "Proteomes", "Dup", prot_name)

	file.copy(clean_prot, clean_dir)
	file.copy(dup_prot, dup_dir)
})
