#!/usr/bin/env Rscript
require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load("circular", "dplyr")

processCOGData		<- function(COG, data = NULL, dataType, penalty, byAge = FALSE, subgroup = NA) {

	if (is.null(data)) {
		stop("Need to define the data object")
	}
	
	if (is.na(penalty)) {
		perCOGAll	<- subset(data[[dataType]]$allPosData, COGcat == COG)
	} else if (!byAge) {
		perCOGAll	<- subset(data[[dataType]][[penalty]]$allPosData, COGcat == COG)
	} else {
		perCOGAll	<- subset(data$lHGT[[penalty]]$allPosData, COGcat == COG & Subgroup == subgroup)
	}

	if (nrow(perCOGAll) == 0) {
		return(NA)
	}

	circStart		<- circular(perCOGAll$relGeneStart * (2 * pi), zero = pi / 2, rotation = "clock", modulo = "2pi")
	perCOGbySubdiv	<- splitBySubdivision(data = perCOGAll, subDivisions = subDivision_list, variable = "COGcat")

	# Combine and extract the data frames containing the per-gene subdivision information. However, this loses the circular data
	subdivColData	<- bind_rows(lapply(perCOGbySubdiv$full_data, function(subdivTable) return(subdivTable$data)))

	# Add the subdivision assignment for each individual gene
	mapCOGontoSubdiv		<- unlist(lapply(perCOGAll$protID, grep, x = subdivColData$protID, fixed = TRUE))
	perCOGAll$SubDivision	<- subdivColData$SubDivision[mapCOGontoSubdiv]
	
	return(list(allData = perCOGAll, CircStart = circStart, perCOGbySubdiv_df = perCOGbySubdiv$tab_df))
}