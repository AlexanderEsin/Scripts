#!/usr/bin/Rscript

## This is currently only applied to long-distance derived HGTs ##

library(ape)
library(phylobase)
library(phytools)
library(stringr)
library(geiger)
library(plyr)
library(gtools)

###########################################################################
## Get the input and output paths ##
input_donor_edge_dir <- "/Users/aesin/Desktop/Mowgli/Long_distance_HGT/Donor_edges/"
input_donor_files <- mixedsort(list.files(input_donor_edge_dir, pattern = "(donor_edges)+"))

output_dir <- "/Users/aesin/Desktop/Geo_analysis/HGT_Donors_basic/Long_distance_HGT/Donor_Taxa/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

###############################
## Get the species tree with node labels ##
species_tree_dir = "/Users/aesin/Desktop/Mowgli/Species_tree/Reconciled_species_tree"

## Convert to a phylo4 object for easy node retrieval ##
setwd(species_tree_dir)
species_tree <- read.tree("outputSpeciesTree.mpr")
species_tree_4 <- phylo4(species_tree)

## Make a df out of all the nodes ##
node_lbs_df <- data.frame(labels(species_tree_4, "all"))
names(node_lbs_df)[1] <- "Labels"

# ## Prepare a df with the edges so we can count ##
# species_edge_df <- data.frame(species_tree$edge)
# names(species_edge_df) <- c("E1", "E2")
# species_edge_df$count <- 0

###############################

setwd(input_donor_edge_dir)
for (input_file in input_donor_files) {
	split_name <- str_split(input_file, "_")[[1]]
	
	penalty <- str_sub(split_name[1], 2)

	if (length(split_name) == 6) {
		min_taxa <- str_sub(split_name[6], 1, -5)
	} else {
		min_taxa <- "0"
	}

	message(paste0("Working on penalty ", penalty, " and a ", min_taxa, " minimum sequence cutoff..."))

	donor_df <- read.table(file = input_file, header = FALSE, sep = "\t")

	## The second edge (edge 2) would be the internal node subtending the tips that are present in that tree ##
	subtended_taxa_list <- vector("list", nrow(donor_df))

	## This works because the labels are listed in node-order in the phylo4 format, meaning that the row-number corresponds to the node number ##
	## NB isntead of taking the lower node, which would truly encompass all subtended taxa - we take the higher node. We do this because Mowgli has a tendency to predict transfers from extinct lineages - if we make the reasonable assumption that the donor group and its sister lineage belong to the same phylum (probably always true), then by taking the ancestral node we increase the chances of capturing some extant taxa ##

	for (i in 1:nrow(donor_df)) {
		## See if it's a tip label first ##
		lower_node <- grep(paste0("_", donor_df[i,2], "$"), node_lbs_df$Labels)
		## If not, assume it's an internal node and search for the specific node number ##
		if (length(lower_node) == 0) {
			lower_node <- grep(paste0("^", donor_df[i,2], "$"), node_lbs_df$Labels)
		}
		## If we can't find a single node, throw an error
		if (length(lower_node) != 1) {
			print("ERROR")
		}

		## Add the group number onto the front of the list ##
		subtended_taxa_list[[i]] <- c(donor_df[i,1], tips(species_tree, lower_node))
		#print(i)
	}

	taxa_list_df <- ldply(subtended_taxa_list, rbind)

	write.table(taxa_list_df, file = paste0(output_dir, "T", penalty, "_refined_donor_taxa_min_", min_taxa, ".txt"), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE, na = "")
}






