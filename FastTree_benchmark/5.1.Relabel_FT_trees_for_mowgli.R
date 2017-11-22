#!/usr/local/bin/Rscript

## This is a shit script that requires work if used again

## Relabel gene tree tips to be used by mowgli. 

## The raxml run produces a tree with tip labels as follows:
## Binomial_name {Class} {Protein_ID.species_ID}
## R automatically collapses the white spaces in the tip labels, so "Binomial_name{Class}{Protein_ID.species_ID}"

## This script removes the class/protein_ID terms, returning just the binomial Binomial_name
## If this binomial name ends in "..._[number]" mowgli would misinterpret this as labelling duplicate taxa, so we replace it with:
## "Binomial_name..._[number]" --> "Binomial_name..._[number]A"

## THEN, for taxa that do appear multiple time, we want to label them as repeats - so we then ADD a "_[number]" onto the end of the tip name


## Packages ##
library(ape)
## Using phyotools midpoint.root() instead of phangorn midpoint()
## midpoint.root() is approx 25x faster
library(phytools)
library(geiger)
library(stringr)

## Set paths
directory		<- "/users/aesin/Desktop/FastTree"
input_trees_dir	<- file.path(directory, "All_trees")
relab_trees_dir	<- file.path(directory, "All_trees_relab")
relab_key_dir	<- file.path(directory, "Tip_keys")

## Make output directories
dir.create(relab_trees_dir, showWarnings = FALSE)
dir.create(relab_key_dir, showWarnings = FALSE)

## A list of the input trees
input_trees <- Sys.glob(file.path(input_trees_dir, "*txt"))

done_counter = 1

for (tree in input_trees) {

	## Extract group number from file name
	tree_file		<- basename(tree)
	group_number	<- str_sub(tree_file, 4, -11)

	## Set the paths for the output files
	output_tree		<- file.path(relab_trees_dir, paste0("Relab_MP_", group_number, ".txt"))
	output_key		<- file.path(relab_key_dir, paste0("Tip_key_", group_number, ".txt"))

	## Read in tree and extract tips
	tree_data		<- read.tree(tree)
	tree_tips		<- tree_data$tip.label

	new_tree		<- tree_data

	trim_tips	<- gsub("\\{.+\\}", "", tree_tips)
	fix_tips	<- str_replace(trim_tips, "(_[0-9]+$)", "\\1A")
	new_tips	<- paste0(fix_tips, "_1")
	tip_df		<- data.frame(Old.Tips = tree_tips, New.Tips = new_tips, stringsAsFactors = FALSE)
	tip_df$Duplicated	<- duplicated(tip_df$New.Tips)

	while (any(tip_df$Duplicated)) {
		current_end	<- as.numeric(str_sub(unique(str_extract(tip_df$New.Tips[which(tip_df$Duplicated == TRUE)], "(_[0-9]+$)")), 2, -1))
		new_end		<- current_end + 1
		tip_df$New.Tips[which(tip_df$Duplicated == TRUE)]	<- str_replace(tip_df$New.Tips[which(tip_df$Duplicated == TRUE)], "(_[0-9]+$)", paste0("_", new_end))
		tip_df$Duplicated	<- duplicated(tip_df$New.Tips)
	}

	## Check that the old tip column in the df is still
	## identically ordered to the tree tip.label
	if (all.equal(new_tree$tip.label, tip_df$Old.Tips) != TRUE) {
		stop("Columns are not equal")
	}

	## Adjust the support values to be out of 100
	new_node_labels	<- floor(as.numeric(new_tree$node.label) * 100)
	new_node_labels[is.na(new_node_labels)] <- "0"

	## Relabel the tips to the new format
	## also relabel the nodes (support values)
	new_tree$tip.label	<- tip_df$New.Tips
	new_tree$node.label	<- new_node_labels
	## Midpoint root the tree
	mp_new_tree 		<- midpoint.root(new_tree)

	## Finally we have to resolve any polychotomies
	mp_new_bi_tree		<- multi2di(mp_new_tree)

	## For the key we just take the old and new tips names as a table
	key_data			<- tip_df[,1:2]

	## Write out the updated tree and key
	write.tree(mp_new_bi_tree, file = output_tree)
	write.table(key_data, file = output_key, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

	message(done_counter)
	done_counter	<- done_counter + 1
}