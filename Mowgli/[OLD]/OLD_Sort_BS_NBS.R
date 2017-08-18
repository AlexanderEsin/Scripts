#!/usr/bin/Rscript

##  THIS HAS NOT BEEN INCORPORATED INTO Relabel_tips_mowgli.R SCRIPT ##

## Take the unsorted output gene tree set and sort it into two directories - those with BS values and those without ##
## The Unsorted_gene_trees directory is simply the renamed Final_BS directory as pulled from the Tree_sorting folder ##

library(ape)

direct <- "/users/aesin/desktop/Mowgli"

unsorted_dir_path <- paste0(direct, "/Unsorted_gene_trees")
sorted_BS_dir_path <- paste0(direct, "/BS_gene_trees")
sorted_NBS_dir_path <- paste0(direct, "/NBS_gene_trees")

dir.create(sorted_BS_dir_path, showWarnings = FALSE)
dir.create(sorted_NBS_dir_path, showWarnings = FALSE)

setwd(unsorted_dir_path)
gene_trees <- dir()

done_counter = 1

for (tree in gene_trees) {

	tree_data <- read.tree(tree)

	BS_output_name <- paste0(sorted_BS_dir_path, "/", tree)
	NBS_output_name <- paste0(sorted_NBS_dir_path, "/", tree)

	if (exists('node.label', where = tree_data) == TRUE) {
		write.tree(tree_data, file = BS_output_name)
	} else {
		write.tree(tree_data, file = NBS_output_name)
	}

	print(done_counter)
	done_counter = done_counter + 1
}

