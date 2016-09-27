#!/usr/bin/Rscript
## Extract the node labels given to anoxy_geobacilli in the species tree ##

library(ape)
library(phangorn)
library(geiger)
library(stringr)

penalty <- 5

parent_dir <- paste0("/Users/aesin/Desktop/Mowgli/Mowgli_outputs/Mow_test_out_t", penalty)
#parent_dir <- paste0("/Users/aesin/Desktop/Mowgli/Test_NBS/Mow_NBS_out_t", penalty)
setwd(parent_dir)
testing_dirs = dir()

counter = 1

for (dir in testing_dirs) {

	setwd(paste0(parent_dir, "/", dir))

	if (file.exists("Anoxy_geo_nodes.tsv")) {
		print(counter)
		counter = counter + 1
		next
	}

	species_tree <- read.tree("outputSpeciesTree.mpr")

	geo_tips <- grep("Geobacillus_", species_tree$tip.label, value = TRUE)
	anoxy_tips <- grep("Anoxybacillus_", species_tree$tip.label, value = TRUE)
	anoxy_geo_tips <- c(geo_tips, anoxy_tips)

	anoxy_geo_prune <- drop.tip(species_tree, setdiff(species_tree$tip.label, anoxy_geo_tips))
	anoxy_geo_nodes <- anoxy_geo_prune$node.label

	for (tip in anoxy_geo_tips) {
		tip_node_label <- str_extract(tip, "[0-9]+$")
		anoxy_geo_nodes <- c(anoxy_geo_nodes, tip_node_label)
	}

	write.table(anoxy_geo_nodes, sep = "\n", file = "Anoxy_geo_nodes.tsv", quote = FALSE, col.names = FALSE, row.names = FALSE)
	print(counter)
	counter = counter + 1
}