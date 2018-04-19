#!/usr/bin/env Rscript

require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load("wesanderson", "dplyr")

# Load master variables and HGT position functions
invisible(sapply(HGTPos.all, source, .GlobalEnv))

# ------------------------------------------------------------------------------------- #
message("Creating the subdivision key...", appendLF = FALSE)

# Compartment ranges (subDivision_list) are defined in the masterVariables.R file

# Colours for the compartments
conPal				<- colorRampPalette(wes_palette("Zissou1"))
subDivison_cols		<- rev(conPal(6))

# ------------------------------------------------------------------------------------- #
# Make a key for the subdivision analysis

# Prepare the data structure
subDivisionKey_df	<- bind_rows(lapply(1:length(subDivision_list), function(half_boundaryIndex) {
	half_boundaries	<- subDivision_list[[half_boundaryIndex]]
	subdiv_name		<- names(subDivision_list)[half_boundaryIndex]
	half_1	<- pi * 2 * half_boundaries
	half_2	<- pi * 2 * (1 - half_boundaries)
	out_df	<- data.frame(xmin = c(half_1[1], half_2[1]), xmax = c(half_1[2], half_2[2]), ymin = rep(0, 2), ymax = rep(1, 2), col = rep(subdiv_name, 2), stringsAsFactors = FALSE)
	return(out_df)
}))

# Factor the colours
subDivisionKey_df$col	<- factor(subDivisionKey_df$col, levels = names(subDivision_list))

message("\rCreating the subdivision key... done")

# ------------------------------------------------------------------------------------- #

# Save the subgroup and position data
message("Saving objects...", appendLF = FALSE)

subDivision_data	<- list(subDivisionKey_df = subDivisionKey_df, subDivison_cols = subDivison_cols)
saveRDS(object = subDivision_data, file = file.path(positionData_path, "subDivisionKeyData.rds"))

message("\rSaving objects... done")
