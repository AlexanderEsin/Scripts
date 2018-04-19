#!/usr/bin/env Rscript

# Load master variables and HGT position functions
invisible(sapply(HGTPos.all, source, .GlobalEnv))

# ------------------------------------------------------------------------------------- #

# Get the branches corresponding to subgroups (for lHGT and sHGT data types)
message("\nProcessing subgroup data...", appendLF = FALSE)

subgroupData		<- getSubgroupData(genome_dir = genome_dir)
subgroupBranch_list	<- subgroupData$subgroupBranch_list

message("\rProcessing subgroup data... done\n")

# ------------------------------------------------------------------------------------- #

# Define the location of the position data
positionInput_path	<- file.path(position_path, "Position_data")

# Define the bandwith for circular density estimation
bandwith 			<- 3000

# Read in position data for different data types
perTypeData	<- lapply(dataTypes, function(dataType) {

	# Position data input directory
	posInput_dir	<- file.path(positionInput_path, paste0(dataType, "_input"))

	# For the "All" gene set, process without penalty
	if (identical(dataType, "All")) return(processInputPositions(dataType = dataType, inputDir = posInput_dir, bandwith = bandwith))

	# For Vertical, lHGT, and sHGT read in the position data per penalty
	perPenaltyData	<- lapply(penalty_list, function(penalty) {
		processed_data	<- processInputPositions(dataType = dataType, penalty = penalty, inputDir = posInput_dir, bandwith = bandwith, subgroupBranches = subgroupBranch_list)
		return(processed_data)
	})

	# Rename the list to penalty values
	names(perPenaltyData)	<- penalty_list

	# Return data for each penalty
	return(perPenaltyData)
})

# Rename list by type of data (e.g "All", "lHGT", etc ..)
names(perTypeData)	<- dataTypes

# ------------------------------------------------------------------------------------- #

# Save the subgroup and position data
message("\nSaving objects...", appendLF = FALSE)

saveRDS(object = subgroupData, file = file.path(positionData_path, "AG_subgroupData.rds"))
saveRDS(object = perTypeData, file = file.path(positionData_path, "AG_perTypeData.rds"))

message("\rSaving objects... done")

# ------------------------------------------------------------------------------------- #
# Done