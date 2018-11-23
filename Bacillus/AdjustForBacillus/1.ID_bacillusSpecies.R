#!/usr/bin/env Rscript

# Any libraries explicitly used in the script
require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load(tidyverse, ape, phangorn)

master		<- file.path("/Users/aesin/Desktop/Bacillus")
tree_dir	<- file.path(master, "Prefilter", "Bacillus_withOut", "Final_trees", "Binomial_labelled")
bac_gen_dir	<- file.path(master, "Bac_genomes", "Genome_gbffs")
bac_list_dir	<- file.path(master, "Bac_genomes", "Genome_lists")
if (!dir.exists(bac_gen_dir)) dir.create(bac_gen_dir, recursive = TRUE)
if (!dir.exists(bac_list_dir)) dir.create(bac_list_dir, recursive = TRUE)

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

edge_tips	<- c("Bacillus_subtilis", "Bacillus_pumilus")
mrca_node	<- getMRCA(mp, edge_tips)

# Extract the clade we want and get the taxids
get_clade	<- extract.clade(mp, node = mrca_node)

# Make ultrametric
ultra_bac	<- chronos(get_clade)

# Cluster this tree
clusters <- ultra_bac %>% 
	dendextend::cutree(k = 25) %>%
	as.matrix(clusters) %>%
	as_tibble(rownames = "ID") %>%
	dplyr::rename(Cluster = 2) %>%
	left_join(taxid_binomd_tbl, by = c("ID" = "Binomial")) %>%
	left_join(acc_ass_taxid_tbl, by = "Taxid") %>%
	left_join(all_genome_data_tbl, by = c("Taxid" = "taxid"))

# Which cluster contains bacillus subtilis str 168
str168	<- clusters %>% filter(str_detect(ID, "str._168"))

# Select all other clusters
non_168_clusts	<- clusters %>%
	filter(!Cluster %in% str168$Cluster)

select_24	<- non_168_clusts %>%
	group_by(Cluster) %>%
	filter(seq_rel_date == max(seq_rel_date))

bacillus_25	<- str168 %>%
	bind_rows(select_24) %>%
	select(ID, Cluster, Taxid, Acc_ass) %>%
	rename(Genome = ID) %>%
	mutate(Genome = str_replace_all(Genome, pattern = "_", replacement = " "))


# Transfer all the relevant gbff files
accession_list <- unlist(lapply(1:nrow(bacillus_25), function(row_index) {
	acc_ass		<- bacillus_25[row_index,] %>% pull(Acc_ass)
	file_name	<- paste0(acc_ass, "_genomic.gbff.gz")
	file_path	<- file.path(master, "Genomes", "Genome_gbffs", file_name)

	file.copy(from = file_path, to = bac_gen_dir, overwrite = TRUE)
	readTop <- read_lines(gzfile(file.path(bac_gen_dir, file_name)), n_max = 10)
	findAcc	<- str_subset(readTop, "ACCESSION")
	procAcc	<- unlist(str_split(findAcc, "\\s+"))[2]

	return(procAcc)
}))

bacillus_25 <- bacillus_25 %>% mutate(Accession = accession_list)

write_tsv(bacillus_25, path = file.path(bac_list_dir, "First_round_select.tsv"))

# ----------------------------------------------------- #
## After checking replichore lengths with dORIc we find that cluster 17 "Bacillus amyloliquefaciens XH7" has a replichore length imbalance
# Remove it and go again with a new genome from cluster 17
non_168_clusts <- non_168_clusts %>% filter(ID != "Bacillus_amyloliquefaciens_XH7")

select_24	<- non_168_clusts %>%
	group_by(Cluster) %>%
	filter(seq_rel_date == max(seq_rel_date))

bacillus_25	<- str168 %>%
	bind_rows(select_24) %>%
	select(ID, Cluster, Taxid, Acc_ass) %>%
	rename(Genome = ID) %>%
	mutate(Genome = str_replace_all(Genome, pattern = "_", replacement = " "))


# Transfer all the relevant gbff files
accession_list <- unlist(lapply(1:nrow(bacillus_25), function(row_index) {
	acc_ass		<- bacillus_25[row_index,] %>% pull(Acc_ass)
	file_name	<- paste0(acc_ass, "_genomic.gbff.gz")
	file_path	<- file.path(master, "Genomes", "Genome_gbffs", file_name)

	file.copy(from = file_path, to = bac_gen_dir, overwrite = TRUE)
	readTop <- read_lines(gzfile(file.path(bac_gen_dir, file_name)), n_max = 10)
	findAcc	<- str_subset(readTop, "ACCESSION")
	procAcc	<- unlist(str_split(findAcc, "\\s+"))[2]

	return(procAcc)
}))

bacillus_25 <- bacillus_25 %>% mutate(Accession = accession_list)

write_tsv(bacillus_25, path = file.path(bac_list_dir, "Second_round_select.tsv"))

# ----------------------------------------------------- #
## This is much better - all replichore lengths are consistently within ~4% of each other
## Now make a list of all the genomes within this Bacillus group that will need to be removed to not interfere with inference
toRemove	<- clusters %>% filter(!Taxid %in% bacillus_25$Taxid)
write_tsv(toRemove, path = file.path(bac_list_dir, "bacillusToRemove.tsv"))

toKeep		<- bacillus_25 %>% mutate(Genome = str_replace_all(Genome, pattern = " ", replacement = "_"))
write_tsv(toKeep, path = file.path(bac_list_dir, "bacillusToKeep.tsv"))

# ----------------------------------------------------- #
# Update new taxid and acc_ass tables for the new dataset

adj_taxid_binomd_tbl	<- taxid_binomd_tbl %>% filter(!Taxid %in% toRemove$Taxid)
adj_acc_ass_taxid_tbl	<- acc_ass_taxid_tbl %>% filter(Taxid %in% adj_taxid_binomd_tbl$Taxid)

write_tsv(adj_taxid_binomd_tbl, path = file.path(bac_list_dir, "bac_TaxidBinomial_table.tsv"))
write_tsv(adj_acc_ass_taxid_tbl, path = file.path(bac_list_dir, "bac_AccAssTaxid_table.tsv"))















