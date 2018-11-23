#!/usr/bin/env Rscript

library(pacman)
p_load(tidyverse, seqinr, magrittr, parallel, ape)

HTa_dir	<- "/Users/aesin/Desktop/HTa_tonio/Round2"

aln_dir	<- file.path(HTa_dir, "Align", "Combined", "COBALT")
tbl_dir	<- file.path(HTa_dir, "Sequences", "Tables")

## ----------------------------------------------- ##

tree_data	<- read.tree(file = file.path(aln_dir, "comb_allOneDom_complete_tree.txt"))
nonBact_tbl	<- read_tsv(file.path(tbl_dir, "allNonBact_table.tsv"), comment = "")

allTips		<- tibble(tipName = tree_data$tip.label) %>%
	left_join(nonBact_tbl %>% select(newID, newTipName), by = c("tipName" = "newID")) %>%
	mutate(replaceTip = case_when(
		!is.na(newTipName) ~ newTipName,
		TRUE ~ tipName))

tree_relab	<- tree_data
tree_relab$tip.label	<- allTips$replaceTip


write.tree(tree_relab, file = file.path(aln_dir, "comb_allOneDom_complete_tree_relab.txt"))