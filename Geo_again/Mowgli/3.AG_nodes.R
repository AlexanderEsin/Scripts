#!/usr/local/bin/Rscript
## Extract the node labels given to anoxy_geobacilli in the species tree ##


if (!require("pacman")) install.packages("pacman")
pacman::p_load("parallel", "ape", "phangorn", "geiger", "stringr", "tools")


extractAGnodes	<- function(specTree, AGtaxids) {
	anoxy_geo_tips	<- grep(paste(paste0(AGtaxids, "_"), collapse = "|"), specTree$tip.label, value = TRUE) 

	anoxy_geo_prune <- drop.tip(specTree, setdiff(specTree$tip.label, anoxy_geo_tips))
	anoxy_geo_nodes <- anoxy_geo_prune$node.label

	tip_nodes		<- str_extract(anoxy_geo_tips, "[0-9]+$")
	anoxy_geo_nodes	<- c(anoxy_geo_nodes, tip_nodes)

	return(anoxy_geo_nodes)
}


################################

mowOut_path <- paste0("/Users/aesin/Desktop/Geo_again/Mowgli/Mowgli_output")
output_dirs	<- dir(mowOut_path, pattern = "Output*", full.names = TRUE)
AGtaxids	<- as.character(read.table(file = "/Users/aesin/Desktop/Geo_again/Genomes/Genome_lists/AG_taxids.txt", sep = "\n")$V1)


################################
## For each penalty, go to the mowgli output dir and within each reconciliation create a list of the anoxy_geo nodes that are necessary downstream ##

numCores	<- detectCores() - 2

AGnodes_extracted	<- lapply(output_dirs, function(output_dir) {

	penalty			<- str_sub(str_extract(basename(output_dir), "_[0-9]"), 2, -1)
	message(paste0("Extracting AG nodes from species tree at penalty: ", penalty, " ..."), appendLF = FALSE)

	groupDirs		<- dir(output_dir, full.names = TRUE)
	
	firstDirTree	<- read.tree(file.path(groupDirs[1], "outputSpeciesTree.mpr"))
	firstDirNodes	<- extractAGnodes(firstDirTree, AGtaxids)
	firstDir_md5	<- md5sum(file.path(groupDirs[1], "outputSpeciesTree.mpr"))

	clust		<- makeCluster(numCores, type = "FORK")

	dirDone		<- parLapply(clust, groupDirs, function(dir) {

		dir_md5		<- md5sum(file.path(dir, "outputSpeciesTree.mpr"))

		if (firstDir_md5 == dir_md5) {
			write.table(firstDirNodes, sep = "\n", file = file.path(dir, "Anoxy_geo_nodes.tsv"), quote = FALSE, col.names = FALSE, row.names = FALSE)
			return("Done")
		}
		dirTree		<- read.tree(file.path(dir, "outputSpeciesTree.mpr"))
		dirNodes	<- extractAGnodes(dirTree, AGtaxids)

		write.table(dirNodes, sep = "\n", file = file.path(dir, "Anoxy_geo_nodes.tsv"), quote = FALSE, col.names = FALSE, row.names = FALSE)
		return("Done")
	})

	message(paste0("\rExtracting AG nodes from species tree at penalty: ", penalty, " ... done"))

	stopCluster(clust)
	return("Done")
})

