#!/usr/local/bin/Rscript

library(ape)
library(phytools)
library(phylobase)
library(phangorn)

NumericTaxidString	<- function(taxid_vector) {
	new_vector	<- paste("tax", taxid_vector, sep = "_")
	return(new_vector)
}

StringTaxidNumeric	<- function(taxid_ST_vector) {
	new_vector	<- str_sub(taxid_ST_vector, 5, -1)
	return(new_vector)
}

direct			<- "/Users/aesin/Desktop/Geo_again"
spec_tree_taxs	<- file.path(direct, "Astral_speciesTree", "Astral_noDup_speciesTree_ultra.txt")
bacill_tree_dir	<- file.path(direct, "Consensus_groups", "Bacillaceae", "Final_trees", "Taxid_labelled")
output_dir		<- file.path(direct, "Astral_speciesTree", "Check_insideGroup")

## Get a list of the Bacillaceae taxids
# Bacillaceae acc_ass and the translation table files
bacillac_acc_ass	<- file.path(direct, "Consensus_groups", "Bacillaceae", "Bacillaceae_acc_ass_list.txt")
anogeo_taxid_file	<- file.path(direct, "Genomes", "Genome_lists", "AG_taxids.txt")
acc_ass_taxid_trans	<- file.path(direct, "Genomes", "Genome_lists", "Acc_ass_taxid_table.tsv")
taxid_binom_trans	<- file.path(direct, "Genomes", "Genome_lists", "Taxid_refinedBinomial_table.tsv")

# Read in the acc_ass and translation files
bac_acc_ass_list	<- scan(bacillac_acc_ass, what = character())
acc_ass_transl		<- read.table(acc_ass_taxid_trans, sep = "\t", header = FALSE, col.names = c("Acc_ass", "Taxid"), colClasses = rep("character", 2))
taxid_bin_transl	<- read.table(taxid_binom_trans, sep = "\t", header = TRUE, colClasses = rep("character", 2))

# Use the acc_ass to get the Bacillaceae taxids
bacillac_taxids		<- acc_ass_transl$Taxid[which(acc_ass_transl$Acc_ass %in% bac_acc_ass_list)]
bacillac_taxids_ST	<- NumericTaxidString(bacillac_taxids)

## Read in Astral species tree and get the subtree containing all Bacillaceae
astral_taxs_data	<- read.tree(spec_tree_taxs)
astral_taxs_MP		<- midpoint(astral_taxs_data)

astral_taxs_MP$tip.label	<- NumericTaxidString(astral_taxs_MP$tip.label)

# ## Now, stupidly, we have to write and re-read this midpoint-rooted tree
# ## Otherwise phylo4 fails to process the tree correctly (WHY???)
# write.tree(astral_taxs_MP, file = file.path(direct, "temp_tree.txt"))
# astral_taxs_MP		<- read.tree(file.path(direct, "temp_tree.txt"))
# file.remove(file.path(direct, "temp_tree.txt"))

## Subset all the taxa from the Astral tree corresponding to the Bacillaceae
## use the most recent common ancestor to see how compactly the Bacillaceae
## are represented in the Astral data.
astral_taxs_MP4		<- phylo4(astral_taxs_MP)
bacillacMRCA		<- MRCA(astral_taxs_MP4, bacillac_taxids_ST)
subset_tree 		<- subset(astral_taxs_MP4, mrca = bacillacMRCA)

## The Bacillaceae are highly associated within the Astral tree - the 215 tips
## are represented fully in a subset tree of 236 tips - i.e. the Bacillaceae are
## not quite monophyletic but the extra taxa pulled out with them appears to also
## be closely related

subset_tips			<- StringTaxidNumeric(as(subset_tree, "phylo")$tip.label)
extra_tips			<- setdiff(subset_tips, bacillac_taxids)
extra_tips_binoms	<- taxid_bin_transl$Binomial[which(taxid_bin_transl$Taxid %in% extra_tips)]

## Let's check the similarity of the monophyletic AnoGeo branching in the Astral
## tree to the branching in the RAxML-reconstructed Bacillaceae tree(s)
anogeo_taxids		<- scan(anogeo_taxid_file, what = character())
anogeo_taxids_ST	<- NumericTaxidString(anogeo_taxids)

tree_evalues		<- list("-100", "-150")
anogeo_rax_trees	<- lapply(tree_evalues, function(trunc_eval) {
	tree_file		<- file.path(bacill_tree_dir, paste0("RAxML_bipartitions.super_tree", trunc_eval, ".txt"))
	tree_data		<- read.tree(tree_file)
	tree_data_MP	<- midpoint(tree_data)
	## Relabel the tips so that phylo4 doesnt get confused
	tree_data_MP$tip.label	<- NumericTaxidString(tree_data_MP$tip.label)
	## Drop the bootstrap values as node labels
	tree_data_4		<- phylo4(tree_data_MP, check.node.labels = "drop")
	anogeoMRCA		<- MRCA(tree_data_4, anogeo_taxids_ST)
	subset_tree		<- subset(tree_data_4, mrca = anogeoMRCA)
	subset_phylo	<- as(subset_tree, "phylo")
	return(subset_phylo)
})
names(anogeo_rax_trees) <- tree_evalues

comp_rax_trees <- cophylo(anogeo_rax_trees$'-100', anogeo_rax_trees$'-150')
plot(comp_rax_trees)

## Extract anogeo from the Astral tree
anogeoMRCA			<- MRCA(astral_taxs_MP4, anogeo_taxids_ST)
anogeo_subset_tree 	<- subset(astral_taxs_MP4, mrca = anogeoMRCA)
anogeo_subset_phylo	<- as(anogeo_subset_tree, "phylo")

comp_astral_100		<- cophylo(anogeo_rax_trees$'-100', anogeo_subset_phylo)
comp_astral_150		<- cophylo(anogeo_rax_trees$'-150', anogeo_subset_phylo)
plot(comp_astral_100)

RF

















