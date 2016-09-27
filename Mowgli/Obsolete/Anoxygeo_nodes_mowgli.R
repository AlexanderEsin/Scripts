#!/usr/bin/Rscript
## Extract the node labels given to anoxy_geobacilli in the species tree ##

library(ape)
library(phangorn)
library(geiger)
library(stringr)

parent_dir = "/Users/aesin/Desktop/Test_dtl_methods/Mowgli/Random_trees_21/Mowgli_output_t4"
setwd(parent_dir)
testing_dirs = dir()

for (dir in testing_dirs) {
	setwd(paste0(parent_dir, "/", dir))

	species_tree <- read.tree("outputSpeciesTree.mpr")

	geo_tips <- grep("Geobacillus_", species_tree$tip.label, value = TRUE)
	anoxy_tips <- grep("Anoxybacillus_", species_tree$tip.label, value = TRUE)
	anoxy_geo_tips <- c(geo_tips, anoxy_tips)

	common_node <- getMRCA(species_tree, anoxy_geo_tips)
	tips(species_tree, common_node)
	anoxy_geo_prune <- drop.tip(species_tree, setdiff(species_tree$tip.label, anoxy_geo_tips))
	anoxy_geo_nodes <- anoxy_geo_prune$node.label

	for (tip in anoxy_geo_tips) {
		tip_node_label <- str_extract(tip, "[0-9]+$")
		anoxy_geo_nodes <- c(anoxy_geo_nodes, tip_node_label)
	}

	write.table(anoxy_geo_nodes, sep = "\n", file = "Anoxy_geo_nodes.tsv", quote = FALSE, col.names = FALSE, row.names = FALSE)
}




# for (node in anoxy_geo_nodes) {
# 	node <- as.character(node)
# 	trans_geo <- which(apply(trans_df, 1, function(x) any(grepl(node,x))))
# 	print(node)
# 	print(trans_geo)
# 	cat("\n")
# }