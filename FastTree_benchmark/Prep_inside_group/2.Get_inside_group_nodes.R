#!/usr/local/bin/Rscript

## 1. Libraries ##
library(ape)
library(phylobase)
library(stringr)

## 2. Paths ##
mowgli_dir	<- "/Users/aesin/Desktop/FastTree"
output_dir	<- file.path(mowgli_dir, "Inside_group")
tree_file	<- file.path(mowgli_dir, "Species_tree", "Reconciled_species_tree", "outputSpeciesTree.mpr")

## 3. Get list of tips. Read in species tree ##
# Read in the list of inside group tips
IG_tip_df	<- read.table(file.path(output_dir, "Inside_group_tips.txt"), header = FALSE, sep = "\n", stringsAsFactors = FALSE)
IG_tip_l	<- split(IG_tip_df, seq(nrow(IG_tip_df)))

# Read in the species tree. We have to use a species tree processed by reconciliation (this assigns the node labels) #
species_tree	<- read.tree(file = tree_file)

## 4. Get mowgli tip names for the inside group tips
IG_tip_renamed_l	<- lapply(IG_tip_l, function(binomial) grep(binomial, species_tree$tip.label, value = TRUE))


## 5. Prune away all branches not associated with IG group - get a list of all the IG nodes ##
IG_prune	<- drop.tip(species_tree, setdiff(species_tree$tip.label, unlist(IG_tip_renamed_l)))
IG_prune_4	<- phylo4(IG_prune)
node_labels	<- as.vector(as.numeric(nodeLabels(IG_prune_4)))

## 6. Write the node labels ##
write.table(node_labels, file = file.path(output_dir, "Inside_group_nodes.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\n")
