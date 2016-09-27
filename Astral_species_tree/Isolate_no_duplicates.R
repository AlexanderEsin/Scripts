#!/usr/bin/Rscript

## This scripts takes a folder full of trees and tests each tree for multiple tips of the same taxon. This paralogy would indicate duplication or HGT without gene replacement or the coalescence of two closely related gene families during the clustering process.
## Trees without any duplicated tip labels are written out to a separate directory.

###########################################################################
## Packages ##
library(ape)
library(phangorn)
library(geiger)
library(stringr)
###########################################################################

## Parent directory ##
direct = "/users/aesin/Desktop/Geo_final_trees"

## The input folder must have labels for taxa that would be identical across all gene trees ##
input_directory = paste0(direct,"/Relabelled")

## Output folder name for the trees not containing duplicates ##
no_duplicates = paste0(direct, "/No_duplicates")
no_duplicates_Geomono = paste0(direct, "/No_duplicates_geo_mono")
dir.create(no_duplicates, showWarnings = FALSE)
dir.create(no_duplicates_Geomono, showWarnings = FALSE)

## Go to input folder and get a list of tree files. Change the glob identifier as necessary ##
setwd(input_directory)
tree_files <- Sys.glob(file.path(getwd(), "Relabelled*"))
# tree_files <- tree_files[1:5]

for (tree in tree_files) {
	
	## Process the name (number) of the tree --> corresponding to the gene family ##
	name = str_extract(tree, "[0-9]+")
	outputname = paste0(no_duplicates, "/Relabelled_", name, ".txt")
	geomono_outputname = paste0(no_duplicates_Geomono, "/Relabelled_", name, ".txt")

	## Read in tree and assign the tips to a variable ##
	tree_data <- read.tree(tree)
	tree_tips <- tree_data$tip.label

	## Make a data frame from all the tips that occur multuple times ##
	n_occur <- data.frame(table(tree_tips))
	multiples <- n_occur[n_occur$Freq > 1, ]

	# print(name)
	# print(nrow(multiples))
	# cat("\n")

	## If the multiples table is empty (there are no duplicated taxa), write the tree to an output folder ##
	if (nrow(multiples) == 0) {
		write.tree(tree_data, file = outputname)
	} else {
		#print(tree)
		next
	}

	## OPTIONAL -- Test whether there Geobacillus is monophyletic and sort those trees that are into a separate folder ##

	## Get all the geobacillus leaf tips ##
	geobac_tips <- grep("Geobacillus_",tree_data$tip.label, value = TRUE)
	if (is.rooted(tree_data) == TRUE) {
		tree_data <- unroot(tree_data)
	}
	## Check whether Geobacillus is monophyletic ##
	if (is.monophyletic(tree_data, geobac_tips)) {
		write.tree(tree_data, file = geomono_outputname)
	}

	## OPTIONAL -- Test whether there are any node labels (BS values in the tree, and sort them into a separate folder) ##

	# if (nrow(multiples) == 0 && exists('node.label', where = tree_data) == TRUE) {
	# 	write.tree(tree_data, file = bs_outputname)
	# } else {
	# 	next
	# }
}