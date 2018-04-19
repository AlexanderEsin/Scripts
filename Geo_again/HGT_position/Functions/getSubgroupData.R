#!/usr/bin/env Rscript
require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load("ape", "phylobase")

getSubgroupData		<- function(genome_dir = NULL) {

	if (is.null(genome_dir)) {
		stop("Need to provide \'Genome\' directory")
	}

	subspeciesGroup_file	<- file.path(genome_dir, "Genome_lists", "AG_subspeciesGroups.txt")
	subspeciesGroup_data	<- read.table(file = subspeciesGroup_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

	# Get those tips that have a fellow tip in the same subgroup. Isolate all the subgroups that have more than one member #
	subspecOnly_data		<- subset(subspeciesGroup_data, duplicated(subspeciesGroup_data$Group) | duplicated(subspeciesGroup_data$Group, fromLast = TRUE))
	uniqueSubGroup_char		<- unique(subspecOnly_data$Group)

	# Read in the taxid <-> mowgliNode translation table (output in Time_analysis.rmd)
	timeAnalysis_dir		<- file.path(master_dir, "HGT_time", "Data")
	taxidMowExtend_df		<- read.table(file = file.path(timeAnalysis_dir, "Taxid2MowTip_table.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
	taxidMowExtend_df$Extension	<- as.character(taxidMowExtend_df$Extension)

	# Read in the consensus AG time tree (output in Time_analysis.rmd)
	AG_conTime_tree			<- read.tree(file = file.path(timeAnalysis_dir, "AG_conTimeTree.tree"))
	AG_conTime_as4			<- phylo4(AG_conTime_tree)

	# Translate the extension-tip tree to binomial-tip tree
	AG_binomTime_tree		<- AG_conTime_tree
	AG_binomTime_tree$tip.label	<- unlist(lapply(AG_binomTime_tree$tip.label, function(tip) {
		newTipLabel_char	<- taxidMowExtend_df$Binomial[which(taxidMowExtend_df$Extension == tip)]
	}))

	# Translate the extension-tip tree to taxid-tip tree
	AG_taxidTime_tree		<- AG_conTime_tree
	AG_taxidTime_tree$tip.label	<- unlist(lapply(AG_taxidTime_tree$tip.label, function(tip) {
		newTipLabel_char	<- taxidMowExtend_df$Taxid[which(taxidMowExtend_df$Extension == tip)]
	}))

	# Work out which branches are subgroup
	subgroup_branches <- vector(mode = "character")
	for (subGroup in uniqueSubGroup_char) {
		taxidsInSubGroup	<- as.vector(subspecOnly_data$Taxid[which(subspecOnly_data$Group == subGroup)])
		taxidsToMowTips		<- taxidMowExtend_df$Extension[which(taxidMowExtend_df$Taxid %in% taxidsInSubGroup)]
		allSubGroupNodes	<- as.vector(descendants(AG_conTime_as4, phylobase::MRCA(AG_conTime_as4, taxidsToMowTips), "all"))

		subgroupEdges_list	<- lapply(allSubGroupNodes, function(node) {
			ancestor	<- ancestor(AG_conTime_as4, node)
			node_lab	<- labels(AG_conTime_as4, type = c("all"))[node]
			ance_lab	<- labels(AG_conTime_as4, type = c("all"))[ancestor]
			branch_lab	<- paste(ance_lab, node_lab, sep = " ")
			return(branch_lab)
		})
		# Add to list of subgroup branches
		subgroup_branches	<- c(subgroup_branches, unlist(subgroupEdges_list))
	}
	return(list(subgroupBranch_list = subgroup_branches, AG_conTime_tree = AG_conTime_tree, AG_taxidTime_tree = AG_taxidTime_tree, AG_binomTime_tree = AG_binomTime_tree))
}
