#!/usr/bin/env Rscript

# Any libraries explicitly used in the script
require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load("dplyr")

# Split data by subdivision (e.g. occurence - for branch) or COG number
splitByZone	<- function(data, zones, variable = NULL) {

	if (identical(nrow(zones), as.integer(8))) {
		column	<- "relGeneStart"
	} else if (identical(nrow(zones), as.integer(5))) {
		column	<- "distToOri"
	} else {
		stop("Don't know which column to use with these zones")
	}

	if (!is.null(variable))	varName	<- unique(unlist(data[[variable]]))

	byZone_data	<- lapply(1:nrow(zones), function(index) {

		zone		<- zones[index,]
		
		# Determine which genes are in each zone
		if (!identical(as.integer(index), nrow(zones))) {
			zone_data			<- data[data[[column]] >= zone$zoneMin & data[[column]] < zone$zoneMax,]
		} else {
			zone_data			<- data[data[[column]] >= zone$zoneMin & data[[column]] <= zone$zoneMax,]
		}

		# Add a column to the raw data being written back out
		zone_data$zoneName	<- if(nrow(zone_data) != 0) zone$zoneName else character(0)
		zone_data$zoneType	<- if(nrow(zone_data) != 0) zone$zoneType else character(0)
		zone_data$zoneCol	<- if(nrow(zone_data) != 0) zone$zoneCol else character(0)

		# Remove the "CircStart" & "CircEnd" columns - we will be combining DFs later on, and binding rows with circular vectors throws warninings
		zone_data				<- zone_data[, !names(zone_data) %in% c("CircStart", "CircEnd")]

		if (is.null(variable)) {
			return(zone_data)
		} else if (identical(variable, "COGcat")) {
			zone_vars			<- unlist(zone_data[[variable]])
			zone_var_tab		<- data.frame(COGcat = varName, zone = zone$zoneName, numObsv = length(zone_vars), stringsAsFactors = FALSE)
		} else {
			zone_vars			<- zone_data[[variable]]
			zone_var_tab		<- data.frame(Variable = varName, zone = zone$zoneName, numObsv = length(zone_vars), stringsAsFactors = FALSE)
		}
		return(list(data = zone_data, cog_tab = zone_var_tab))
	})

	if (is.null(variable)) {
		byZoneCombined_data	<- bind_rows(byZone_data)
		return(byZoneCombined_data)
	} else {
		# Extract the tables and bind them into one dataframe
		zoneCombine_tables		<- lapply(byZone_data, function(zone) return(zone$cog_tab))
		zoneCombine_tables_df	<- bind_rows(zoneCombine_tables)
		zoneCombine_tables_df	<- zoneCombine_tables_df %>% group_by(zone) %>% summarise(numObsv = sum(numObsv), COGcat = unique(zoneCombine_tables_df$COGcat))
		# Any missing values to 0
		zoneCombine_tables_df[is.na(zoneCombine_tables_df)]	<- 0
		# Add proportion column
		zoneCombine_tables_df$Proportion	<- zoneCombine_tables_df$numObsv / sum(zoneCombine_tables_df$numObsv)
		# Full data
		byZoneCombined_data	<- bind_rows(lapply(byZone_data, function(zone) return(zone$data)))

		return(list(full_data = byZoneCombined_data, tab_df = zoneCombine_tables_df))
	}
}