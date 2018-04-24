#!/usr/bin/env Rscript

# Any libraries explicitly used in the script
require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load("dplyr")

# Split data by subdivision (e.g. occurence - for branch) or COG number
splitBySubdivision	<- function(data, subDivisions, variable = NA) {

	data$Mod_dist	<- data$relGeneStart
	data$Mod_dist[which(data$Mod_dist > 0.5)]	<- 1 - data$Mod_dist[which(data$Mod_dist > 0.5)]
	varName			<- unique(unlist(data[[variable]]))

	bysubdiv_data	<- lapply(1:length(subDivisions), function(index) {

		subdiv		<- subDivisions[[index]]
		subdivName	<- names(subDivisions[index])


		if (!identical(index, length(subDivisions))) {
			subdiv_data			<- data[which(data$Mod_dist >= subdiv[1] & data$Mod_dist < subdiv[2]),]	
		} else {
			subdiv_data			<- data[which(data$Mod_dist >= subdiv[1] & data$Mod_dist <= subdiv[2]),]
		}

		# Add a column to the raw data being written back out
		subdiv_data$SubDivision	<- if(nrow(subdiv_data) != 0) subdivName else character(0)

		# Remove the "CircStart" & "CircEnd" columns - we will be combining DFs later on, and binding rows with circular vectors throws warninings
		subdiv_data				<- subset(subdiv_data, select = c(-CircStart, -CircEnd))

		if (is.na(variable)) {
			return(data.frame(subDiv = subdivName, numObsv = nrow(subdiv_data), stringsAsFactors = FALSE))
		} else if (identical(variable, "COGcat")) {
			subdiv_vars			<- unlist(subdiv_data[[variable]])
			subdiv_var_tab		<- data.frame(COGcat = varName, subDiv = subdivName, numObsv = length(subdiv_vars), stringsAsFactors = FALSE)
		} else {
			subdiv_vars			<- subdiv_data[[variable]]
			subdiv_var_tab		<- data.frame(Variable = varName, subDiv = subdivName, numObsv = length(subdiv_vars), stringsAsFactors = FALSE)
		}
		return(list(data = subdiv_data, cog_tab = subdiv_var_tab))
	})

	if (is.na(variable)) {
		subdivCombine_df	<- bind_rows(bysubdiv_data)
		return(subdivCombine_df)
	} else {
		# Extract the tables and bind them into one dataframe
		subdivCombine_tabs		<- lapply(bysubdiv_data, function(subdiv) return(subdiv$cog_tab))
		subdivCombine_df		<- bind_rows(subdivCombine_tabs)
		# Any missing values to 0
		subdivCombine_df[is.na(subdivCombine_df)]	<- 0
		# Add proportion column
		subdivCombine_df$Proportion	<- subdivCombine_df$numObsv / sum(subdivCombine_df$numObsv)

		return(list(full_data = bysubdiv_data, tab_df = subdivCombine_df))
	}
}
