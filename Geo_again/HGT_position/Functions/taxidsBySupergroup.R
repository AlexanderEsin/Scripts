#!/usr/bin/env Rscript
require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load("parallel")


taxidsBySupergroup	<- function(supergroup, speciesGroupings = NULL, nodes_data = NULL, names_data = NULL) {

	# Check all necessary data is provided
	if (is.null(speciesGroupings)) stop("Provide subgroupings data by species: output in Genome_lists/Species_groupings.tsv")
	if (is.null(nodes_data) || is.null(names_data)) stop("Provide both nodes and names taxdmp data")

	# What are we working on
	message(paste0("Identifying taxids belonging to ", supergroup, "..."))

	# Get the taxid and rank of the supergroup
	supergroupTaxid	<- names_data$taxid[which(names_data$name == supergroup & names_data$nameClass == "scientific name")]
	supergroupRank	<- nodes_data$nodeRank[which(nodes_data$taxid == supergroupTaxid)]

	# Unique parental nodes in species groupings
	uniqueParentSpecies		<- unique(speciesGroupings$parentSpecies)

	# For each parent species, find whether it belongs to this supergroup
	superGroupTaxids_list	<- mclapply(uniqueParentSpecies, function(taxid) {
		
		childTaxid	<- taxid

		while (TRUE) {
			parentEntry	<- nodes_data[which(nodes_data$taxid == childTaxid),]

			if (identical(parentEntry$parentTaxid, supergroupTaxid)) {
				return(taxid)
			} else if (identical(parentEntry$nodeRank, supergroupRank) || identical(parentEntry$parentTaxid, "1")) {
				break
			} else {
				childTaxid	<- parentEntry$parentTaxid
			}
		}
		return(NA)
	}, mc.cores = 20)

	# Remove all empty elements (not part of supergroup) and return
	superGroupTaxids	<- unlist(superGroupTaxids_list[!is.na(superGroupTaxids_list)])
	return(superGroupTaxids)
}