#!/usr/bin/env Rscript

# # Load master variables and HGT position functions
# invisible(sapply(HGTPos.all, source, .GlobalEnv))

# require(pacman, warn.conflicts = FALSE, quietly = TRUE)
# p_load("dplyr", "ggplot2", "wesanderson")

# # ------------------------------------------------------------------------------------- #`
# # Read in data

# message("\nReading in data...", appendLF = FALSE)

# # General position data
# perTypeData			<- readRDS(file.path(positionData_path, "AG_perTypeData.rds"))

# # Subgroup data
# subgroupData		<- readRDS(file.path(positionData_path, "AG_subgroupData.rds"))

# # Subdivision key
# subDivisionKey_file	<- file.path(positionData_path, "subDivisionKeyData.rds")
# if(!file.exists(subDivisionKey_file)) {
# 	stop("The subdivision data file \"", subDivisionKey_file, "\" is required.\nRun the subDivisionPrep.R script!")
# } else {
# 	subDivisionKey_data	<- readRDS(file.path(positionData_path, "subDivisionKeyData.rds"))
# }

# message("\rReading in data... done\n")

# # ------------------------------------------------------------------------------------- #
# # Output path for figures
# GIanalysisFig_path	<- file.path(figureOutput_path, "genomicIslands")
# if (!dir.exists(GIanalysisFig_path)) dir.create(GIanalysisFig_path)

# # ------------------------------------------------------------------------------------- #
# # Set global quartz options
# quartz.options(canvas = "#333233", bg = "#333233")

# # ------------------------------------------------------------------------------------- #
# # Read in GI data
# giProcess_path		<- file.path(master_path, "HGT_position", "GI_data")
# giBoundary_file		<- file.path(giProcess_path, "perGenome_relGIBoundaries.tsv")
# giBoundary_data		<- read.table(giBoundary_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# GI_presence_list	<- lapply(1:nrow(giBoundary_data), function(index) {

# 	GI_row	<- giBoundary_data[index,]
# 	if (is.na(GI_row$GI_relStart)) {
# 		return(NA)
# 	}

# 	GI_present	<- seq(from = round(GI_row$GI_relStart, 4), to = round(GI_row$GI_relEnd, 4), by = 0.0001)
# 	out_df	<- data.frame(Species = rep(GI_row$Binomial, length(GI_present)), GI_location = GI_present, stringsAsFactors = FALSE)
# 	return(out_df)
# })

# # Remove all the genomes with no GI data
# GI_presence_list	<- GI_presence_list[!is.na(GI_presence_list)]
# GI_presence_df		<- bind_rows(GI_presence_list)

# # Number of unique (GI present) species
# numSpeciesGI		<- length(unique(GI_presence_df$Species))

# # Overall plot
# GI_crossGenomeDens_plot	<- ggplot(data = GI_presence_df, aes(x = GI_location)) +
# 	scale_x_continuous(limits = c(0, 1), name = "Relative Genome Position") +
# 	geom_histogram(aes(y = ..ncount.. / 4, fill = Species), bins = 1000) +
# 	geom_density(adjust = 0.25, color = "#D9D9D9", show.legend = FALSE) +
# 	scale_fill_manual(values = alpha(wes_palette("Darjeeling1", n = numSpeciesGI, type = "continuous"), 0.8)) +
# 	darkTheme


# quartz(width = 20, height = 10, type = "pdf", file = file.path(GIanalysisFig_path, "GI_crossGenomePositions.pdf"))
# print(GI_crossGenomeDens_plot)
# invisible(dev.off())



# # Remove the two species with genome rearrangements
# rearrangedSpecies	<- c("Anoxybacillus amylolyticus", "Geobacillus sp. Y412MC61")
# numSpeciesNoRearr	<- numSpeciesGI - length(rearrangedSpecies)

# # Plot without rearranged genomes
# GI_noRearrGenomes_plot	<- ggplot(data = subset(GI_presence_df, !Species %in% rearrangedSpecies), aes(x = GI_location)) +
# 	scale_x_continuous(limits = c(0, 1), name = "Relative Genome Position") +
# 	scale_y_continuous(name = "Density of GI regions") +
# 	geom_histogram(aes(y = ..ncount.. / 4, fill = Species), bins = 1000) +
# 	geom_density(adjust = 0.25, color = "#D9D9D9", show.legend = FALSE) +
# 	scale_fill_manual(values = alpha(wes_palette("Darjeeling1", n = numSpeciesNoRearr, type = "continuous"), 0.8)) +
# 	darkTheme


# # ------------------------------------------------------------------------------------- #
# # Process data to separate GI and non-GI data points

# # 19 species have annotated GIs
# speciesWithGI	<- unique(GI_presence_df$Species)

# byType_withGI		<- lapply(list("lHGT", "sHGT", "Ver"), function(dataType) {

# 	message(paste0("Data Type = ", dataType))

# 	byPenalty_withGI	<- lapply(penalty_list, function(penalty) {

# 		message(paste0("\tWorking on penalty = ", penalty))

# 		bySpeciesWithGI_list	<- mclapply(speciesWithGI, function(species) {

# 			message(species)

# 			bySpeciesPos_data	<- subset(perTypeData[[dataType]][[penalty]]$allPosData, binomial == species)
# 			bySpeciesGI_data	<- subset(giBoundary_data, Binomial == species)

# 			withinGI	<- unlist(lapply(1:nrow(bySpeciesPos_data), function(geneIndex) {

# 				geneData		<- bySpeciesPos_data[geneIndex,]
# 				relGeneStart	<- geneData$relGeneStart
				
# 				withinAnyGI		<- unlist(lapply(1:nrow(bySpeciesGI_data), function(GI_index) {
# 					giEntry		<- bySpeciesGI_data[GI_index,]
# 					withinGI	<- ifelse(relGeneStart >= giEntry$GI_relStart & relGeneStart <=giEntry$GI_relEnd, TRUE, FALSE)
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
# names(byType_withGI)	<- c("lHGT", "sHGT", "Ver")


# ------------------------------------------------------------------------------------- #
# Compare the disribution of all lHGT genes, lHGTs in GIs and lHGTS outside of GIs

lHGT_subsets	<- c("All", "In_GI", "Out_GI")
lHGT_penalty	<- "4"
bandwith		<- 3000

# All genes for species with annotated GIs
allPosData_withGI		<- subset(perTypeData$All$allPosData, binomial %in% speciesWithGI)
allPosData_withGI_dens	<- density.circular(allPosData_withGI$CircStart, kernel = "vonmises", bw = bandwith)

# Set up plot
quartz(width = 21, height = 8)
par(mfrow = c(1, 3))
par(mar = c(0, 0, 0, 0))

invisible(lapply(lHGT_subsets, function(subset) {

	# Are we looking for genes within GIs?
	getInGI	<- ifelse(identical(subset, "In_GI"), TRUE, FALSE)

	# Subset data depending on whether we want genes in GIs
	if (identical(subset, "All")) {
		subset_data	<- byType_withGI$lHGT[[lHGT_penalty]]
	} else {
		subset_data	<- subset(byType_withGI$lHGT[[lHGT_penalty]], In_GI == getInGI)
	}

	# Calculate circular density
	subset_dens	<- density.circular(subset_data$CircStart, kernel = "vonmises", bw = bandwith)

	# Produce circular plot
	circularDensityPlot(dataDensityA = subset_dens,
		bgDensity = allPosData_withGI_dens,
		shrink = 0.7,
		tcl.offset = 0.8,
		titleCex = 1.8,
		titleName = paste0(subset, "\nGenes = ", nrow(subset_data)))
}))

outputFileName	<- file.path(GIanalysisFig_path, "GIvsNon-GI_lHGT_density.pdf")
invisible(dev.copy2pdf(file = outputFileName))
invisible(dev.off())


# ------------------------------------------------------------------------------------- #
# Compare the Outside GI lHGT gene distribution across multiple penalties


# Set up plot
quartz(width = 12, height = 12)
par(mfrow = c(2, 2))
par(mar = c(0, 0, 0, 0))

invisible(lapply(penalty_list, function(penalty) {

	# Subset data by penalty
	subset_data	<- subset(byType_withGI$lHGT[[penalty]], In_GI == FALSE)

	# Calculate circular density
	subset_dens	<- density.circular(subset_data$CircStart, kernel = "vonmises", bw = bandwith)

	# Produce circular plot
	circularDensityPlot(dataDensityA = subset_dens,
		bgDensity = allPosData_withGI_dens,
		shrink = 0.7,
		tcl.offset = 0.8,
		titleCex = 1.8,
		titleName = paste0("Penalty = ", penalty, "\nGenes = ", nrow(subset_data)))
}))

outputFileName	<- file.path(GIanalysisFig_path, "crossPenaltyNon-GI_density.pdf")
invisible(dev.copy2pdf(file = outputFileName))
invisible(dev.off())



# ------------------------------------------------------------------------------------- #
# Ran a test with the downloaded GI data for GKau (with locus IDs from Island Viewer)
# They use an earlier assembly revision, so 4 proteins we have they are lacking. Thus, 93/97 proteins overlap

gkauGIViewer_data	<- read.table(file.path(giProcess_path, "testGKau.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# lHGTinGI_crossCheck	<- testlHGT_data[which(testlHGT_data$locusTag %in% gkauGIViewer_data$Locus),]
# lHGTinGI_myMeth[which(!lHGTinGI_myMeth$locusTag %in% lHGTinGI_crossCheck$locusTag),]































