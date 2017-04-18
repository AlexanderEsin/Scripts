#!/usr/local/bin/Rscript
### In this script: check that Anoxy/Geobacillus with vertical descent are monophyletic in the gene tree ###

## Libraries ##
library(ape)
library(plyr)
library(gtools)

## Directories / files for the input list and trees ##
input_file <- "/Users/aesin/Desktop/Geo_analysis/New_geobacillus_consensus/Vertical_family_prot_ids.tsv"
tree_directory <- "/Users/aesin/Desktop/Geo_analysis/Geo_v_all/2.0/Trees/All_trees/"

## Read in table ##
family_prots_tbl <- read.csv(file = input_file, sep = "\t", header = FALSE, colClasses = "character")

## Function to find the tip in the gene tree based on protein ID ##
FindTipInTree <- function(prot_id, tree_tips) {
	tip_in_tree <- grep(prot_id, tree_tips, value = TRUE)
}

## Function to check whether all the protein IDs are monophyletic in the tree ##
CheckMonophyletic <- function(family) {

	# Split input into: gene family number, all Anoxy/Geobacillus protein IDs #
	fam_num		<- family[1]
	fam_prots	<- family[2:length(family)]

	# Construct tree file name #
	tree_name <- paste0(tree_directory, fam_num, "_tree.txt")
	# file.exists(tree_name); # only run as a test once

	# Read in tree, and isolate tip labels #
	gene_tree <- read.tree(file = tree_name)
	gene_tree_tips <- gene_tree$tip.label

	# Find the correct tip in the tree based on the protein ID #
	anoxygeo_tree_tips <- as.character(lapply(fam_prots, FindTipInTree, tree_tips = gene_tree_tips))

	# Make sure all tips found - no errors #
	if (length(anoxygeo_tree_tips) != 22) {
		message(paste0("Tree ", fam_num, " does not have 22 tips."))
	}

	# Check if monophyletic #
	monophyly_state <- is.monophyletic(phy = gene_tree, tips = anoxygeo_tree_tips, reroot = FALSE)

	# Make the group/fam number as character (otherwise output looks bad) #
	group_num <- as.character(fam_num)

	return(c(group_num, monophyly_state))
}

## Produces slightly stupid output but nicely formatted ##
raw_results <- ddply(family_prots_tbl, 1, CheckMonophyletic)

## Remove all NAs and concatenate all the TRUE / FALSE into a single column. Remove intermediate columns ##
names(raw_results) <- c("Gene_family", "Monophyly")

## Sort correctly ##
results <- raw_results[mixedorder(raw_results$Gene_family),]
rownames(results) <- NULL

## Write output as just gene families that have all 22 species, and are monophyletic ##
write.table(results[,1], file = "/Users/aesin/Desktop/Geo_analysis/New_geobacillus_consensus/Vertical_monophyletic_fams.tsv", sep = "\n", quote = FALSE, row.names = FALSE, col.names = FALSE)
