#!/usr/local/bin/Rscript
## Extract the node labels given to anoxy_geobacilli in the species tree ##

require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load(tidyverse, fs, ape, phangorn, parallel, geiger, tools)
 
extractCoreNodes	<- function(specTree, coreTaxids) {

	taxidRelab		<- str_subset(specTree$tip.label, paste(paste0("^", coreTaxids, "_"), collapse = "|"))

	coreMRCA		<- getMRCA(specTree, taxidRelab)
	corePrune_tree	<- extract.clade(specTree, coreMRCA)

	intNodes		<- corePrune_tree$node.label
	tipNodes		<- str_extract(taxidRelab, "[0-9]+$")

	allCoreNodes	<- tibble(coreNode = c(intNodes, tipNodes))

	return(allCoreNodes)
}


################################

master		<- "/Users/aesin/Desktop/Staph"
mowOut_path <- file.path(master, "Mowgli", "Mowgli_output")
mow_dirList	<- dir(mowOut_path, pattern = "Output*", full.names = TRUE)

coreKeep_file	<- file.path(master, "Core_genomes", "Genome_lists", "coreToKeep.tsv")
coreKeep_tbl	<- read_tsv(coreKeep_file)


################################
## For each penalty, go to the mowgli output dir and within each reconciliation create a list of the anoxy_geo nodes that are necessary downstream ##

coreNodes_extracted	<- lapply(mow_dirList, function(outDir) {

	penalty			<- str_extract(basename(outDir), "(?<=_)[0-9]+")
	message(paste0("Extracting Core nodes from species tree at penalty: ", penalty, " ..."), appendLF = FALSE)

	groupDirs		<- dir(outDir, full.names = TRUE)
	firstDirTree_f	<- file.path(groupDirs[1], "outputSpeciesTree.mpr")
	
	firstDirTree	<- read.tree(firstDirTree_f)
	firstDirNodes	<- extractCoreNodes(firstDirTree, coreKeep_tbl$Taxid)
	firstDir_md5	<- md5sum(firstDirTree_f)

	invisible(mclapply(groupDirs, function(groupDir) {

		dirTree_f	<- file.path(groupDir, "outputSpeciesTree.mpr")
		dir_md5		<- md5sum(dirTree_f)

		if (firstDir_md5 == dir_md5) {
			write_tsv(firstDirNodes, path = file.path(groupDir, "coreNodes.tsv"), col_names = FALSE)
			return(1)
		}

		if (identical(file.size(dirTree_f), 0)) {
			dir_delete(groupDir)
			return(0)
		}

		dirTree		<- read.tree(dirTree_f)
		dirNodes	<- extractCoreNodes(dirTree, coreKeep_tbl$Taxid)

		write_tsv(dirNodes, path = file.path(groupDir, "coreNodes.tsv"), col_names = FALSE)
		return(1)
	}, mc.cores = 10))

	message(paste0("\rExtracting Core nodes from species tree at penalty: ", penalty, " ... done"))
	return("Done")
})

