#!/usr/bin/Rscript

## This script identfies a common leaf set between a set of trees and a given reference tree and prunes the gene trees down to the common leaf set so that they it can be directly compared to the referene tree
## The output is pruned trees that may contain duplicate taxa names as a result of paralogy

###########################################################################
## Packages ##
library(ape)
library(phangorn)
library(geiger)
library(stringr)

###########################################################################

## Parent directory ##
direct = "/users/aesin/Desktop/Test_import_group"

## Input and output directories ##
input_directory = paste0(direct, "/Relabelled")
reference_directory = paste0(direct)
pruned_directory = paste0(direct, "/Pruned")
create.dir(pruned, showWarnings = FALSE)

## Get tip names from the reference tree ##
setwd(reference_directory)
reference_tree = Sys.glob(file.path(getwd(), "Reference*"))
reference_tree_data <- read.tree(reference_tree)
reference_tree_tips <- reference_tree_data$tip.label

## Go to input folder and get a list of tree files. Change the glob identifier as necessary ##
setwd(input_directory)
tree_files <- Sys.glob(file.path(getwd(), "Relabelled*"))

for (tree in tree_files) {

	## Process the name (number) of the tree --> corresponding to the gene family ##
	name = str_extract(tree, "[0-9]+")
	outputname = paste0(pruned_directory, "/Pruned_", name, ".nwk")

	## Read in tree and assign the tips to a variable ##
	tree_data <- read.tree(tree)
	tree_tips <- tree_data$tip.label

	## Get a common set of tips that are shared between the reference tree and the gene tree ##
	shared_tips = intersect(reference_tree_tips, tree_tips)

	## If the number of shared tips is fewer than 4, cannot produce a pruned tree ##
	if (length(shared_tips) < 4) {
		next
	}

	## Prune the gene tree to contain only those tips that are shared with the reference tree ##
	pruned_tree <- drop.tip(tree_data, setdiff(tree_tips, shared_tips))

	## OPTIONAL: If there are any duplicated tips in the resultant pruned tree, notify by calling infomration to STDOUT ## 
	if (anyDuplicated(pruned_tree$tip.label) > 0) {
		cat("\nDuplicated tips:", name, "Inside tips:", length(shared_tips))
		next
	}
	
	## Write the pruned tree to the pruned_directory ##
	write.tree(pruned_tree, file = outputname)
}