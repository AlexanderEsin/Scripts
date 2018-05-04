#!/usr/bin/env Rscript

# Load master variables and HGT position functions
invisible(sapply(HGTPos.all, source, .GlobalEnv))

# ------------------------------------------------------------------------------------- #
# Subgroup information

# Get the branches corresponding to subgroups (for lHGT and sHGT data types)
message("\nProcessing subgroup data...", appendLF = FALSE)

subgroupData		<- getSubgroupData(genome_dir = genome_path, databasePath = allProtDB_path)
subgroupBranch_list	<- subgroupData$subgroupBranch_list

message("\rProcessing subgroup data... done\n")

# ------------------------------------------------------------------------------------- #
# Position data for AG

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
# Data structure for COG analysis
message("\nProcessing COG data... ")

# List of all the unique COG categorties found in the entire AG dataset
uniqueCOGs		<- unique(unlist(perTypeData$All$allPosData$COGcat))

perTypeCOG_data	<- lapply(dataTypes_withAge, function(dataType) {

	message(paste0("\tProcessing COG data for \'", dataType, "\'... "), appendLF = FALSE)

	byAge		<- FALSE
	subgroup	<- NA 

	if (identical(dataType, "All")) {
		penalty		<- NA
	} else if (identical(dataType, "Ver")) {
		penalty		<- verPenalty
	} else {
		penalty		<- hgtPenalty
		if (identical(dataType, "Old") || identical(dataType, "Recent")) {
			byAge		<- TRUE
			subgroup	<- ifelse(identical(dataType, "Old"), FALSE, TRUE)
		}
	}

	# For each COG, extract the relevant data - see the "processCOGData" function for exact outputs
	perCOGData_list			<- lapply(uniqueCOGs, processCOGData, data = perTypeData, dataType = dataType, penalty = penalty, byAge = byAge, subgroup = subgroup)
	names(perCOGData_list)	<- uniqueCOGs
	
	# Remove any COGs that don't have any observations
	perCOGData_list			<- perCOGData_list[!is.na(perCOGData_list)]

	# Produce the circular plots
	perCOGAvPos_Sum_Plot	<- COGdistribCircular_plot(perCOGData_list = perCOGData_list, dataType = dataType)
	perCOGAvPos_Sum			<- perCOGAvPos_Sum_Plot$circSum
	perCOGAvPos_Plot		<- perCOGAvPos_Sum_Plot$plot

	# Produce the by-Compartment dataframe
	byCOGbySubdiv_df		<- bind_rows(lapply(perCOGData_list, function(COG) return(COG$perCOGbySubdiv_df)))
	byCOGbySubdiv_df$Set	<- dataType

	message(paste0("\r\tProcessing COG data for \'", dataType, "\'... done"))

	# Return all the raw datam as well as COGs per compartment, the circular summary of gene positions, and the circular weighted mean position plot
	return(list(perCOGdata = perCOGData_list, bySubdivision = byCOGbySubdiv_df, circSummary = perCOGAvPos_Sum, plot = perCOGAvPos_Plot))
})
names(perTypeCOG_data)	<- dataTypes_withAge


# ------------------------------------------------------------------------------------- #

# Save the subgroup and position data
message("\nSaving objects...", appendLF = FALSE)

saveRDS(object = subgroupData, file = file.path(positionData_path, "AG_subgroupData.rds"))
saveRDS(object = perTypeData, file = file.path(positionData_path, "AG_perTypeData.rds"))
saveRDS(object = perTypeCOG_data, file = file.path(positionData_path, "AG_perTypeCOGData.rds"))

message("\rSaving objects... done")

# ------------------------------------------------------------------------------------- #
# Done