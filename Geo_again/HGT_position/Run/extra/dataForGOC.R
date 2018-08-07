#!/usr/bin/env Rscript

# Load master variables and HGT position functions
invisible(sapply(HGTPos.all, source, .GlobalEnv))

# Any libraries explicitly used in the script
require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load("tidyverse")

# ------------------------------------------------------------------------------------- #

# Read in all position data
perTypeData			<- readRDS(file.path(positionData_path, "AG_perTypeData.rds"))

# ------------------------------------------------------------------------------------- #

# Extract vertical data at penalty = 3
allVertical_data	<- perTypeData$Ver$'3'$allPosData

oneToOneOrths_only	<- allVertical_data %>%
	# Group by orthologous groups
	group_by(orthGroup) %>%
	# If there are 23 binomials in group and all 23 are unique, it's a one-to-one orth
	summarise(one2one_Orth = case_when(
			length(unique(binomial)) == 23 & length(binomial) == 23 ~ TRUE,
			TRUE ~ FALSE
		)
	) %>%
	# Select only one-to-one orthologs
	filter(one2one_Orth == TRUE) %>%
	# Vector of ortholog groups numbers
	pull(orthGroup) %>%
	# Filter initial data frame to contain only one-to-one groups
	{filter(allVertical_data, orthGroup %in% .)} %>%
	# Remove unwanted columns
	select(-c(CircStart, CircEnd))

# Split the dataframe into list - per binomial
splitBySpecies	<- split(oneToOneOrths_only, oneToOneOrths_only$binomial)

# ------------------------------------------------------------------------------------- #

# List structure can be converted back into single dataframe (not run):
# require(dplyr)
# x	<- bind_rows(splitBySpecies)

# Save as RDS object
saveRDS(splitBySpecies, file = file.path(positionData_path, "OneToOne_VerticalOrthologs_23Genomes.rds"))