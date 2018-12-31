#!/usr/bin/env Rscript

# Any libraries explicitly used in the script
require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load(tidyverse, fs, ape, phangorn)

master			<- file.path("/Users/aesin/Desktop/Staph")
grp_list_dir	<- file.path(master, "Core_genomes", "Genome_lists")
tree_dir		<- file.path(master, "Prefilter", "Staphylococcus", "Final_trees", "Binomial_labelled")

outTree_dir		<- file.path(master, "Core_tree")
dir_create(outTree_dir)

toKeep_tbl	<- read_tsv(file.path(grp_list_dir, "coreToKeep.tsv"))


# Process tree
binom_file	<- file.path(tree_dir, "FT_binomial_SupTree_-50.txt")
binom_data	<- read.tree(binom_file)
binom_mp	<- midpoint(binom_data)

prune_tree	<- keep.tip(binom_mp, tip = toKeep_tbl$Genome)
write.tree(prune_tree, file.path(outTree_dir, "coreOnly_tree.txt"))
