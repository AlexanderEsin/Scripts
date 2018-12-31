#!/usr/bin/env Rscript

# Load master variables and HGT position functions
invisible(sapply(HGTPos.bac, source, .GlobalEnv))

# ------------------------------------------------------------------------------------- #
# Subgroup information

# Get the branches corresponding to subgroups (for lHGT and sHGT data types)
# message("\nProcessing subgroup data...", appendLF = FALSE)

# subgroupData		<- getSubgroupData(genome_dir = genome_path, databasePath = allProtDB_path)
# subgroupBranch_list	<- subgroupData$subgroupBranch_list

# message("\rProcessing subgroup data... done\n")

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
	if (identical(dataType, "All")) return(processInputPositions(dataType = dataType, inputDir = posInput_dir, bandwith = bandwith, removeSpecies = outlierTaxid))

	# For Vertical, lHGT, and sHGT read in the position data per penalty
	perPenaltyData	<- lapply(penalty_list, function(penalty) {
		processed_data	<- processInputPositions(dataType = dataType, penalty = penalty, inputDir = posInput_dir, bandwith = bandwith, removeSpecies = outlierTaxid)
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

# Zone list
zoneBoundary_file	<- file.path(positionData_path, "AG_zoneBoundaries.rds")
if (!file.exists(zoneBoundary_file)) stop("Run the 3.zoneIdentification.R script first")

zoneBoundaryList	<- readRDS(zoneBoundary_file)

# List of all the unique COG categorties found in the entire AG dataset
uniqueCOGs		<- unique(unlist(perTypeData$All$allPosData$COGcat))

perTypeCOG_data	<- lapply(dataTypes, function(dataType) {

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
	perCOGData_list			<- lapply(uniqueCOGs, processCOGData, data = perTypeData, dataType = dataType, penalty = penalty, byAge = byAge, subgroup = subgroup, zone_data = zoneBoundaryList$fullRange)
	names(perCOGData_list)	<- uniqueCOGs
	
	# Remove any COGs that don't have any observations
	perCOGData_list			<- perCOGData_list[!is.na(perCOGData_list)]

	# Produce the by-Compartment dataframe
	byCOGbyZone_df		<- bind_rows(lapply(perCOGData_list, function(COG) return(COG$perCOGbyZone_df)))
	byCOGbyZone_df$Set	<- dataType

	message(paste0("\r\tProcessing COG data for \'", dataType, "\'... done"))

	# Return all the raw datam as well as COGs per compartment, the circular summary of gene positions, and the circular weighted mean position plot
	return(list(perCOGdata = perCOGData_list, byZone = byCOGbyZone_df))
})
names(perTypeCOG_data)	<- dataTypes


# # ------------------------------------------------------------------------------------- #
# # Process GI data
# message("\nProcessing GI data... ")

# # Read in GI boundary file
# giBoundary_file		<- file.path(giProcess_path, "perGenome_relGIBoundaries.tsv")
# giBoundary_data		<- read.table(giBoundary_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# # Partition the genome into segments covered by GI
# GI_presence_list	<- lapply(1:nrow(giBoundary_data), function(index) {

# 	GI_row	<- giBoundary_data[index,]
# 	if (is.na(GI_row$GI_relStart)) {
# 		return(NA)
# 	}

# 	GI_present	<- seq(from = round(GI_row$GI_relStart, 4), to = round(GI_row$GI_relEnd, 4), by = 0.0001)
# 	out_df	<- data.frame(
# 		taxid = rep(GI_row$Taxid, length(GI_present)),
# 		binomial = rep(GI_row$Binomial, length(GI_present)),
# 		GI_location = GI_present,
# 		stringsAsFactors = FALSE)
# 	return(out_df)
# })

# # Remove all the genomes with no GI data
# GI_presence_list	<- GI_presence_list[!is.na(GI_presence_list)]
# GI_presence_df		<- bind_rows(GI_presence_list)

# # Remove the outlier species
# GI_positions_data		<- subset(GI_presence_df, !taxid %in% outlierTaxid)

# # Prepare a data object for the GI positions to be written out
# GI_positions_out		<- list(rawBoundary_data = giBoundary_data, GI_positions_data = GI_positions_data)

# # Number of unique (GI present) species
# speciesWithGI	<- unique(GI_positions_data$binomial)
# numSpeciesGI	<- length(speciesWithGI)

# # For each gene, does it lie within the GI boundaries for that genome?
# byType_withGI_data		<- lapply(dataTypes, function(dataType) {

# 	message(paste0("Data Type = ", dataType))

# 	if (identical(dataType, "All")) {

# 		bySpeciesWithGI_list	<- mclapply(speciesWithGI, function(species) {			

# 			bySpeciesPos_data	<- subset(perTypeData$All$allPosData, binomial == species)
# 			bySpeciesGI_data	<- subset(giBoundary_data, Binomial == species)

# 			withinGI	<- unlist(lapply(1:nrow(bySpeciesPos_data), function(geneIndex) {

# 				geneData		<- bySpeciesPos_data[geneIndex,]
# 				relGeneStart	<- geneData$relGeneStart
				
# 				withinAnyGI		<- unlist(lapply(1:nrow(bySpeciesGI_data), function(GI_index) {
# 					giEntry		<- bySpeciesGI_data[GI_index,]
# 					withinGI	<- ifelse(relGeneStart >= giEntry$GI_relStart & relGeneStart <= giEntry$GI_relEnd, TRUE, FALSE)
# 					return(withinGI)
# 				}))

# 				if(any(withinAnyGI)) return(TRUE) else return(FALSE)
# 			}))

# 			# Remove the circStart and circEnd columns as we'll break the circular class by bind_rows
# 			bySpeciesPos_data	<- subset(bySpeciesPos_data, select = -c(CircStart, CircEnd))

# 			bySpeciesWithGI_data	<- cbind(bySpeciesPos_data, In_GI = withinGI, stringsAsFactors = FALSE)
# 			return(bySpeciesWithGI_data)
# 		}, mc.cores = 14)

# 		# Dataframe including only species with annotated GIs - each entry has a TRUE/FALSE for whether the gene is within a GI
# 		bySpeciesWithGI_df		<- bind_rows(bySpeciesWithGI_list)

# 		# Circularise start and end positions
# 		bySpeciesWithGI_df$CircStart	<- circular(bySpeciesWithGI_df$relGeneStart * (2 * pi), zero = pi / 2, rotation = "clock", modulo = "2pi")
# 		bySpeciesWithGI_df$CircEnd		<- circular(bySpeciesWithGI_df$relGeneEnd * (2 * pi), zero = pi / 2, rotation = "clock", modulo = "2pi")

# 		return(bySpeciesWithGI_df)
# 	}
	
# 	byPenalty_withGI	<- lapply(penalty_list, function(penalty) {

# 		message(paste0("\tWorking on penalty = ", penalty))

# 		bySpeciesWithGI_list	<- mclapply(speciesWithGI, function(species) {

			

# 			bySpeciesPos_data	<- subset(perTypeData[[dataType]][[penalty]]$allPosData, binomial == species)
# 			bySpeciesGI_data	<- subset(giBoundary_data, Binomial == species)

# 			withinGI	<- unlist(lapply(1:nrow(bySpeciesPos_data), function(geneIndex) {

# 				geneData		<- bySpeciesPos_data[geneIndex,]
# 				relGeneStart	<- geneData$relGeneStart
				
# 				withinAnyGI		<- unlist(lapply(1:nrow(bySpeciesGI_data), function(GI_index) {
# 					giEntry		<- bySpeciesGI_data[GI_index,]
# 					withinGI	<- ifelse(relGeneStart >= giEntry$GI_relStart & relGeneStart <= giEntry$GI_relEnd, TRUE, FALSE)
# 					return(withinGI)
# 				}))

# 				if(any(withinAnyGI)) return(TRUE) else return(FALSE)
# 			}))

# 			# Remove the circStart and circEnd columns as we'll break the circular class by bind_rows
# 			bySpeciesPos_data	<- subset(bySpeciesPos_data, select = -c(CircStart, CircEnd))

# 			bySpeciesWithGI_data	<- cbind(bySpeciesPos_data, In_GI = withinGI, stringsAsFactors = FALSE)
# 			return(bySpeciesWithGI_data)
# 		}, mc.cores = 14)

# 		# Dataframe including only species with annotated GIs - each entry has a TRUE/FALSE for whether the gene is within a GI
# 		bySpeciesWithGI_df		<- bind_rows(bySpeciesWithGI_list)

# 		# Circularise start and end positions
# 		bySpeciesWithGI_df$CircStart	<- circular(bySpeciesWithGI_df$relGeneStart * (2 * pi), zero = pi / 2, rotation = "clock", modulo = "2pi")
# 		bySpeciesWithGI_df$CircEnd		<- circular(bySpeciesWithGI_df$relGeneEnd * (2 * pi), zero = pi / 2, rotation = "clock", modulo = "2pi")

# 		return(bySpeciesWithGI_df)
# 	})
# 	names(byPenalty_withGI)	<- penalty_list
# 	return(byPenalty_withGI)
# })
# names(byType_withGI_data)	<- dataTypes




# ------------------------------------------------------------------------------------- #
# Save the subgroup and position data
message("\nSaving objects...", appendLF = FALSE)

# saveRDS(object = subgroupData, file = file.path(positionData_path, "AG_subgroupData.rds"))
saveRDS(object = perTypeData, file = file.path(positionData_path, "AG_perTypeData.rds"))
saveRDS(object = perTypeCOG_data, file = file.path(positionData_path, "AG_perTypeCOGData.rds"))
# saveRDS(object = byType_withGI_data, file = file.path(positionData_path, "AG_perTypeGIData.rds"))
# saveRDS(object = GI_positions_out, file = file.path(positionData_path, "AG_GI_positions.rds"))

message("\rSaving objects... done")
