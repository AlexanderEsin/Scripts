#!/usr/local/bin/Rscript
## Extract the node labels given to anoxy_geobacilli in the species tree ##

library(ape)
library(phangorn)
library(geiger)
library(stringr)
library(tools)


################################
## Script variables ##
penalty_list	<- list(6)
parent_dir		<- "/Users/aesin/Desktop/FastTree/FastTree_outputs"

################################
## For each penalty, go to the mowgli output dir and within each reconciliation create a list of the anoxy_geo nodes that are necessary downstream ##

for (penalty in penalty_list) {

	penalty_dir		<- file.path(parent_dir, paste0("Output_", penalty))
	message(paste0("Now working on penalty = ", penalty))


	testing_dirs = dir(penalty_dir, full.names = TRUE)

	counter = 1

	for (dir in testing_dirs) {

		if (counter == 1) {
			md5_model	<- md5sum(file.path(dir, "outputSpeciesTree.mpr"))
		} else {
			md5_tree	<- md5sum(file.path(dir, "outputSpeciesTree.mpr"))
			if (md5_model == md5_tree) {

				write.table(model_nodes, sep = "\n", file = file.path(dir, "Anoxy_geo_nodes.tsv"), quote = FALSE, col.names = FALSE, row.names = FALSE)
				
				print(counter)
				counter = counter + 1
				next
			}
		}

		species_tree	<- read.tree(file.path(dir, "outputSpeciesTree.mpr"))

		# if (file.exists("Anoxy_geo_nodes.tsv")) {
		# 	print(counter)
		# 	counter = counter + 1
		# 	next
		# }

		match_names		<- c("Geobacillus_", "Anoxybacillus_")
		anoxy_geo_tips	<- grep(paste(match_names, collapse = "|"), species_tree$tip.label, value = TRUE) 

		anoxy_geo_prune <- drop.tip(species_tree, setdiff(species_tree$tip.label, anoxy_geo_tips))
		anoxy_geo_nodes <- anoxy_geo_prune$node.label

		tip_nodes		<- str_extract(anoxy_geo_tips, "[0-9]+$")
		anoxy_geo_nodes	<- c(anoxy_geo_nodes, tip_nodes)

		if (counter == 1) {
			model_nodes	<- anoxy_geo_nodes
		}

		write.table(anoxy_geo_nodes, sep = "\n", file = file.path(dir, "Anoxy_geo_nodes.tsv"), quote = FALSE, col.names = FALSE, row.names = FALSE)
		print(counter)
		counter = counter + 1
	}

}