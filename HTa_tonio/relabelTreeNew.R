#!/usr/bin/env Rscript

library(pacman)
p_load(tidyverse, ape, parallel)

## ----------------------------------------------- ##

HTa_dir		<- "/Users/aesin/Desktop/HTa_tonio/Round2"
txd_dir		<- file.path(HTa_dir, "taxdmp")
tree_dir	<- file.path(HTa_dir, "RAxML")
tbl_dir		<- file.path(HTa_dir, "Sequences", "Tables")

# Taxon files
nodes_data	<- as_tibble(read_rds(file.path(txd_dir, "nodes.rds")))
names_data	<- as_tibble(read_rds(file.path(txd_dir, "names.rds")))

# bact ID tbl
bact_ID_tbl	<- as_tibble(read_tsv(file.path(tbl_dir, "byFamSubsetBact_select.tsv"))) %>%
	mutate(newName = str_replace_all(newName, " ", "_"))
nonBact_ID_tbl	<- as_tibble(read_tsv(file.path(tbl_dir, "allNonBact_table.tsv"))) %>%
	select(newID, newTipName)
# raxml tree
rawRax_tree	<- read.tree(file.path(tree_dir, "RAxML_bipartitions.complete_aln_relab_raxml_tree.txt"))

## ----------------------------------------------- ##

tipLabels	<- tibble(rawTree_tips = rawRax_tree$tip.label) %>%
	mutate(trimTree_tips = str_replace_all(rawTree_tips, "_Full$", ""))

# Extract just the bacterial tips
bactOnly	<- tipLabels %>%
	filter(trimTree_tips %in% bact_ID_tbl$newName) %>%
	left_join(bact_ID_tbl, by = c("trimTree_tips" = "newName"))

# Identify the phylum for the bacteria
phylumData_df	<- mclapply(1:nrow(bactOnly), function(index) {

	row			<- bactOnly[index,]
	taxidNum	<- row$familyTaxid

	while (TRUE) {
		if (identical(taxidNum, 1)) break

		entry	<- nodes_data %>% filter(taxid == taxidNum)
		rank	<- entry$nodeRank
		if (identical(rank, "phylum")) break
		taxidNum	<- entry$parentTaxid
	}

	if (identical(taxidNum, 1)) {
		message(paste0("Index = ", index, ". PHYLUM NOT FOUND!"))
		return(tibble(phylumTaxid = NA, phylumName = NA))
	}

	phylumTaxid	<- entry$taxid
	phylumName	<- names_data %>% filter(taxid == phylumTaxid & nameClass == "scientific name") %>% pull(name)
	out_tbl		<- tibble(phylumTaxid = phylumTaxid, phylumName = phylumName)

	message(paste0("Index = ", index, ". ", phylumName))

	return(out_tbl)
}, mc.cores = 10)
phylumData_df	<- bind_rows(phylumData_df)
bact_wPhylum	<- bind_cols(bactOnly, phylumData_df)

# Where phylum cannot be found, use "NA [familyName]"
phylumFix		<- bact_wPhylum %>%
	mutate(phylumFix = case_when(
		!is.na(phylumName) ~ phylumName,
		TRUE ~ familyName)) %>%
	group_by(phylumFix) %>%
	mutate(phylumTip = paste0(phylumFix, "_", seq(1:length(phylumFix)))) %>%
	ungroup()

# Select subset for joining
bact_wPhy_prune	<- phylumFix %>% select(rawTree_tips, phylumTip)

newTipLabels	<- tipLabels %>%
	left_join(bact_wPhy_prune, by = "rawTree_tips") %>%
	left_join(nonBact_ID_tbl, by = c("rawTree_tips" = "newID")) %>%
	mutate(newTreeTips = case_when(
		is.na(phylumTip) ~ newTipName,
		!is.na(phylumTip) ~ phylumTip))

write_tsv(phylumFix, path = file.path(tbl_dir, "byPhylum_subset_bactTable.tsv"))


tipRelab_tree	<- rawRax_tree
tipRelab_tree$tip.label	<- newTipLabels$newTreeTips


write.tree(tipRelab_tree, file.path(tree_dir, "relabelTip_raxmlTree.txt"))



