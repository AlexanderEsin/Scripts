#!/usr/bin/env Rscript
require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load("stringr", "circular")

processInputPositions	<- function(dataType = "All", penalty = NA, inputDir, bandwith, subgroupBranches = NA) {

	# Track progress
	message(paste0("\tProcessing \'", dataType, "\' at penalty: ", penalty, " ..."), appendLF = FALSE)

	# Set the name of the input file depending on data type
	if (!identical(dataType, "All")) {
		if (is.na(penalty))  stop("Provide a penalty unless dataType = \'All\'")
		fileID	<- paste0("T", penalty, "_", dataType)
	} else {
		fileID	<- dataType
	}

	# Find the file, check it exists and read it in
	genePos_file	<- file.path(inputDir, paste0(fileID, "_positionData.tsv"))
	if (!file.exists(genePos_file)) stop(paste0("Cannot find file ", genePos_file))
	genePos_data	<- read.table(file = genePos_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

	# Circularise the relative start and end positions
	genePos_data$CircStart	<- circular(genePos_data$relGeneStart * (2 * pi), zero = pi / 2, rotation = "clock", modulo = "2pi")
	genePos_data$CircEnd	<- circular(genePos_data$relGeneEnd * (2 * pi), zero = pi / 2, rotation = "clock", modulo = "2pi")

	# Calculate circular density (on the start position)
	genePosDens		<- density.circular(genePos_data$CircStart, kernel = "vonmises", bw = bandwith)

	# For the HGT sets, determine which transfers are into subgroup (vs group) branches
	if (identical(dataType, "lHGT") || identical(dataType, "sHGT")) {
		if (is.na(penalty))  stop("Need subgroup branch data for HGT dataTypes")
		genePos_data$Subgroup	<- FALSE
		genePos_data$Subgroup[which(genePos_data$recepEdge %in% subgroupBranches)]	<- TRUE
	}

	# Process to COG column - make each entry a list of n COGs
	genePos_data$COGcat[is.na(genePos_data$COGcat)]	<- "-"
	genePos_data$COGcat	<- str_split(genePos_data$COGcat, "\\|")

	# Track progress and return
	message(paste0("\rProcessing \'", dataType, "\' at penalty: ", penalty, " ... done"))

	return(list(Penalty = penalty, allPosData = genePos_data, circDensity = genePosDens))
}

