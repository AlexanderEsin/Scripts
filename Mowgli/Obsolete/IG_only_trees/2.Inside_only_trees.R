#!/usr/bin/Rscript
## Get a list of trees that contain only Inside Group taxa ##

library(ape)
library(phangorn)
library(geiger)
library(stringr)
library(gtools)

###########################################################################

## Get the original insde group tips ##
setwd("/Users/aesin/Desktop/Mowgli")
in_tips_file <- "Inside_group_tips.txt"
inside_group_tips <- scan(in_tips_file, character(0), quiet = TRUE)

## Read the original species tree ##
setwd("/Users/aesin/Desktop/Mowgli/Species_tree")
orig_species_tree <- read.tree("Ultrametric_species_tree_tipfix.txt")

correct_in_tips = character()
for (in_tip in inside_group_tips) {
	correct_tip <- grep(in_tip, orig_species_tree$tip.label, value = TRUE)
	correct_in_tips = c(correct_in_tips, correct_tip)
}

###########################################################################

setwd("/Users/aesin/Desktop/Mowgli/Final_BS")
all_bs_trees <- dir()
inside_only_trees = character()

done_counter = 1

for (tree in all_bs_trees) {

	gene_tree <- read.tree(tree)

	for (n in 1:length(gene_tree$tip.label)) {

	  ## Remove the bracketed terms ##
	  old_label <- gene_tree$tip.label[n]
	  new_label <- gsub("\\{.+\\}", "", old_label)

	  ## Attempt to fix any trailing "_[number]" labels ##
	  num_fix_label <- str_replace(new_label, "(_[0-9]+$)", "\\1A")

	  ## Check whether the num_fix_label is the same as the original trimmed label - i.e. was there a trailing number that was modified? If so, make the fixed label the new output ##
	  if (identical(new_label, num_fix_label) == FALSE) {
	    new_label <- num_fix_label
	  }

	  gene_tree$tip.label[n] <- new_label
	}

	total_tips <- gene_tree$tip.label
	x <- setdiff(total_tips, correct_in_tips)

	## If all the gene tree tips correspond to insde group tips, then the difference will be 0 ##
	if (length(x) == 0) {
		tree_number = str_extract(tree, "[0-9]+")
		inside_only_trees = c(inside_only_trees, tree_number)
	}
	print(done_counter)
	done_counter = done_counter + 1
}

inside_only_trees <- mixedsort(inside_only_trees)

setwd("/Users/aesin/Desktop/Mowgli")
write.table(as.matrix(inside_only_trees), file = "Trees_in_group_only.tsv", row.names = FALSE, col.names = FALSE, quote = FALSE)
