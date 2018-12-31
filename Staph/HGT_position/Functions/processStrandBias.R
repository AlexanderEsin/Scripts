#!/usr/bin/env Rscript
require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load("dplyr", "zoo", "parallel")

processStrandBias	<- function(dataType, data = NULL, verPenalty = "3", hgtPenalty = "4", binNumber = 100) {

	if (is.null(data)) {
		stop("Need to define the data object")
	}

	message(paste0("Calculating strand bias profiles for: ", dataType))

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

	byTypebySpeciesStrandBias	<- mclapply(binomial_list, function(binomialName) {

		byBinStrandBias		<- lapply(seq(1:binNumber), function(genome_bin) {

			bin_start		<- (genome_bin - 1) / binNumber
			bin_end			<- genome_bin / binNumber

			if (is.na(penalty)) {
				perBinStrand	<- subset(data$All$allPosData, binomial == binomialName & relGeneStart > bin_start & relGeneStart <= bin_end, select = relStrand, drop = TRUE)
			} else if (!byAge) {
				perBinStrand	<- subset(data[[dataType]][[penalty]]$allPosData, binomial == binomialName & relGeneStart > bin_start & relGeneStart <= bin_end, select = relStrand, drop = TRUE)
			} else {
				perBinStrand	<- subset(data$lHGT[[penalty]]$allPosData, binomial == binomialName & Subgroup == subgroup & relGeneStart > bin_start & relGeneStart <= bin_end, select = relStrand, drop = TRUE)
			}

			if (length(perBinStrand) == 0) {
				return(data.frame(BinIndex = bin_end, StrandBias = NA, stringsAsFactors = FALSE))
			}

			numSameStrand	<- length(which(perBinStrand == "same"))
			numDiffStrand	<- length(which(perBinStrand == "diff"))

			sameOvDiff		<- (numSameStrand - numDiffStrand) / sum(c(numSameStrand, numDiffStrand))
			return(data.frame(BinIndex = bin_end, StrandBias = sameOvDiff, stringsAsFactors = FALSE))
		})
		byBinStrandBias_df	<- bind_rows(byBinStrandBias)
		byBinStrandBias_df$Species	<- binomialName

		## Add a rolling mean value to identify outlier patterns in the "All" dataType
		if (identical(dataType, "All"))  {
			byBinStrandBias_df$rollMean	<- as.vector(rollapply(zoo(byBinStrandBias_df$StrandBias), width = 20, mean, by = 1, fill = NA))
		}
		return(byBinStrandBias_df)
	}, mc.cores = 20)
	byTypebySpeciesStrandBias		<- bind_rows(byTypebySpeciesStrandBias)
	byTypebySpeciesStrandBias$Type	<- dataType
	return(byTypebySpeciesStrandBias)
}

