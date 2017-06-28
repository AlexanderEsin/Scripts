library(ape)
library(phylobase)
library(phytools)
library(stringr)
library(geiger)
library(plyr)

###########################################################################

penalty = 4
min_taxa = 0

input_donor_edge_dir <- "/Users/aesin/Desktop/Mowgli/External_HGT/New_parser/Donor_edges/"

if (min_taxa > 4) {
	file_name_add <- paste0("_min_", min_taxa)
} else {
	file_name_add <- ""
}
###########################################################################

species_tree_dir = "/Users/aesin/Desktop/Mowgli/Species_tree/Sample_mowgli_species_tree"
setwd(species_tree_dir)

species_tree <- read.tree("outputSpeciesTree.mpr")
species_tree_4 <- phylo4(species_tree)

node_lbs_df <- data.frame(labels(species_tree_4, "all"))
names(node_lbs_df)[1] <- "Labels"

## Prepare a data frame with the edges so we can count ##
species_edge_df <- data.frame(species_tree$edge)
names(species_edge_df) <- c("E1", "E2")
species_edge_df$count <- 0

donor_df <- read.table(file = paste0(input_donor_edge_dir, "T", penalty, "_refined_donor_edges", file_name_add, ".txt"), header = FALSE, sep = "\t")

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
}

taxa_list_df <- ldply(subtended_taxa_list, rbind)
write.table(taxa_list_df, file = paste0(input_donor_edge_dir, "T", penalty, "_refined_donor_taxa_min_", min_taxa), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE, na = "")