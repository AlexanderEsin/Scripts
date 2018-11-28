#!/usr/bin/env Rscript

# Any libraries explicitly used in the script
require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load(tidyverse, fs, ape, phangorn)

master		<- file.path("/Users/aesin/Desktop/Staph")
tree_dir	<- file.path(master, "Prefilter", "Staphylococcus", "Final_trees", "Binomial_labelled")
grp_gen_dir	<- file.path(master, "Staph_genomes", "Genome_gbffs")
grp_list_dir	<- file.path(master, "Staph_genomes", "Genome_lists")
dir_create(c(grp_gen_dir, grp_gen_dir))

# Translate taxid / binomial
taxid_binomd_tbl	<- read_tsv(file.path(master, "Genomes", "Genome_lists", "Taxid_refinedBinomial_table.tsv"))
acc_ass_taxid_tbl	<- read_tsv(file.path(master, "Genomes", "Genome_lists", "Acc_ass_taxid_table.tsv"), col_names = c("Acc_ass", "Taxid"))
all_genome_data_tbl	<- read_tsv(file.path(master, "Genomes", "Genome_lists", "All_complete_genomes.tsv")) %>%
	select(taxid, refseq_category, version_status, assembly_level, seq_rel_date)


# ----------------------------------------------------- #
# Process tree
binom_file	<- file.path(tree_dir, "FT_binomial_SupTree_-50.txt")
binom_data	<- read.tree(binom_file)

# Midpoint root
mp	<- midpoint(binom_data)
# Confirm the root as expected
# plot(mp)

edge_tips	<- c("Staphylococcus_aureus", "Staphylococcus_stepanovicii")
mrca_node	<- getMRCA(mp, edge_tips)

# Extract the clade we want and get the taxids
get_clade	<- extract.clade(mp, node = mrca_node)

# Make ultrametric
ultra_clade	<- chronos(get_clade)

# Cluster this tree
clusters <- ultra_clade %>% 
	dendextend::cutree(k = 25) %>%
	as.matrix(clusters) %>%
	as_tibble(rownames = "ID") %>%
	dplyr::rename(Cluster = 2) %>%
	left_join(taxid_binomd_tbl, by = c("ID" = "Binomial")) %>%
	left_join(acc_ass_taxid_tbl, by = "Taxid") %>%
	left_join(all_genome_data_tbl, by = c("Taxid" = "taxid"))

# Which cluster contains Staphylococcus_aureus_subsp._aureus_NCTC_8325 
refStrain	<- clusters %>% filter(str_detect(ID, "NCTC"))

# Select all other clusters
non_ref_clusts	<- clusters %>%
	filter(!Cluster %in% refStrain$Cluster)

select_24	<- non_ref_clusts %>%
	group_by(Cluster) %>%
	filter(seq_rel_date == max(seq_rel_date))

staph_25	<- refStrain %>%
	bind_rows(select_24) %>%
	select(ID, Cluster, Taxid, Acc_ass) %>%
	dplyr::rename(Genome = ID) %>%
	mutate(Genome = str_replace_all(Genome, pattern = "_", replacement = " "))


# Transfer all the relevant gbff files
accession_list <- unlist(lapply(1:nrow(staph_25), function(row_index) {
	acc_ass		<- staph_25[row_index,] %>% pull(Acc_ass)
	file_name	<- paste0(acc_ass, "_genomic.gbff.gz")
	file_path	<- file.path(master, "Genomes", "Genome_gbffs", file_name)

	file.copy(from = file_path, to = grp_gen_dir, overwrite = TRUE)
	readTop <- read_lines(gzfile(file.path(grp_gen_dir, file_name)), n_max = 10)
	findAcc	<- str_subset(readTop, "ACCESSION")
	procAcc	<- unlist(str_split(findAcc, "\\s+"))[2]

	return(procAcc)
}))

staph_25 <- staph_25 %>% mutate(Accession = accession_list)

write_tsv(staph_25, path = file.path(grp_list_dir, "First_round_select.tsv"))


# ----------------------------------------------------- #
## Now make a list of all the genomes within this group that will need to be removed to not interfere with inference
toRemove	<- clusters %>% filter(!Taxid %in% staph_25$Taxid)
write_tsv(toRemove, path = file.path(grp_list_dir, "grpToRemove.tsv"))

toKeep		<- staph_25 %>% mutate(Genome = str_replace_all(Genome, pattern = " ", replacement = "_"))
write_tsv(toKeep, path = file.path(grp_list_dir, "grpToKeep.tsv"))

# ----------------------------------------------------- #
# Update new taxid and acc_ass tables for the new dataset

adj_taxid_binomd_tbl	<- taxid_binomd_tbl %>% filter(!Taxid %in% toRemove$Taxid)
adj_acc_ass_taxid_tbl	<- acc_ass_taxid_tbl %>% filter(Taxid %in% adj_taxid_binomd_tbl$Taxid)

write_tsv(adj_taxid_binomd_tbl, path = file.path(grp_list_dir, "grp_TaxidBinomial_table.tsv"))
write_tsv(adj_acc_ass_taxid_tbl, path = file.path(grp_list_dir, "grp_AccAssTaxid_table.tsv"))















