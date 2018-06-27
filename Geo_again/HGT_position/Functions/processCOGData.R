#!/usr/bin/env Rscript
require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load("circular", "dplyr")

processCOGData		<- function(COG, data = NULL, dataType, penalty, byAge = FALSE, subgroup = NA, zone_data = NULL) {

	if (is.null(data)) {
		stop("Need to define the data object")
	} else if (is.null(zone_data)) {
		stop("Need to provide zone data")
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
	perCOGbyZone	<- splitByZone(data = perCOGAll, zones = zone_data, variable = "COGcat")

	return(list(allData = perCOGbyZone$full_data, CircStart = circStart, perCOGbyZone_df = perCOGbyZone$tab_df))
}