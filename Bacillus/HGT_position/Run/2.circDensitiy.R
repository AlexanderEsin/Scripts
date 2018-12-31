#!/usr/bin/env Rscript

# Load master variables and HGT position functions
invisible(sapply(HGTPos.bac, source, .GlobalEnv))

require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load(tidyverse, wesanderson, gridExtra)

# ------------------------------------------------------------------------------------- #
# Read in data
message("\nReading in data...", appendLF = FALSE)

perTypeData				<- readRDS(file.path(positionData_path, "AG_perTypeData.rds"))

message("\rReading in data... done\n")

# ------------------------------------------------------------------------------------- #  
# Output path for figures
circDensityFig_path	<- file.path(figureOutput_path, "circDensity")
if (!dir.exists(circDensityFig_path)) dir.create(circDensityFig_path)

# ------------------------------------------------------------------------------------- #
# Set global quartz options
quartz.options(canvas = "white", bg = "white")



# ------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------- #


# if (!is.null(outlierTaxid)) {
	# Pick bandwith
	snglPlot_bw	<- 3000

	# Plot the density for a single penalty - manuscript
	HGT_data	<- perTypeData$lHGT$'4'$allPosData
	Ver_data	<- perTypeData$Ver$'3'$allPosData
	all_data	<- perTypeData$All$allPosData

	# Calculate densities
	HGT_density	<- density.circular(HGT_data$CircStart, kernel = "vonmises", bw = snglPlot_bw)
	Ver_density	<- density.circular(Ver_data$CircStart, kernel = "vonmises", bw = snglPlot_bw)
	All_density	<- density.circular(all_data$CircStart, kernel = "vonmises", bw = snglPlot_bw)

	# Number of genes for each density plot


	# Plotting this is tricky: straight-to-pdf gives split output. Plot to quartz device and read off device to PDF file
	quartz(width = 12, height = 12)

	par(mar = c(0,0,0,0))

	circularDensityPlot(
		dataDensityA = HGT_density,
		dataDensityB = Ver_density,
		bgDensity = All_density,
		uin = 3,
		shrink = 0.5,
		tcl.offset = 0.8,
		bg = "white",
		axisCol = axisCol,
		enrichUpColor = dataTypeCols$HGT,
		enrichDownColor = wes_palette("Darjeeling1")[1],
		densBLineColor = dataTypeCols$Ver)


	outputFileName	<- file.path(circDensityFig_path, "singlePenalty_HGTDensity_bw3000_plot.pdf")
	invisible(dev.copy2pdf(file = outputFileName))
	invisible(dev.off())
# }




# ------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------- #
# Vertical gene density across different penalty values
message("Plotting Vertical gene density across penalties...", appendLF = FALSE)

quartz(width = 14, height = 14)
par(mfrow = c(2, 2))
par(mar = c(0, 0, 0, 0))
custom_bw	<- 8000
invisible(lapply(penalty_list, function(penalty) {

	# Density values needed
	if (is.null(custom_bw)) {
		verDensity	<- perTypeData$Ver[[penalty]]$circDensity
		allDensity	<- perTypeData$All$circDensity
	} else {
		message(paste0("Using custom bandwidth to re-estimate circular density. Bandwith = ", custom_bw))
		ver_circPos	<- perTypeData$Ver[[penalty]]$allPosData$CircStart
		all_circPos	<- perTypeData$All$allPosData$CircStart

		verDensity	<- density.circular(ver_circPos, kernel = "vonmises", bw = custom_bw)
		allDensity	<- density.circular(all_circPos, kernel = "vonmises", bw = custom_bw)
	}
	
	# Number of genes for each density plot
	numGenes	<- length(verDensity$data)

	# Produce plot
	position_plot	<- circularDensityPlot(dataDensityA = verDensity,
		bgDensity = allDensity,
		shrink = 0.9,
		tcl.offset = 0.8,
		titleName = paste0("Vertical Enrichment at Penalty = ", penalty, "\nGenes = ", numGenes))
	replayPlot(position_plot)
}))
outputFileName	<- file.path(circDensityFig_path, "VerbyPenalty.pdf")
invisible(dev.copy2pdf(file = outputFileName))
invisible(dev.off())

message("\rPlotting Vertical gene density across penalties... done")

# ------------------------------------------------------------------------------------- #
# lHGT density across different penalty values
message("Plotting lHGT gene density across penalties...", appendLF = FALSE)

quartz(width = 12, height = 12)
par(mfrow = c(2, 2))
par(mar = c(0, 0, 0, 0))
invisible(lapply(penalty_list, function(penalty) {

	# Density values needed
	lHGTDensity	<- perTypeData$lHGT[[penalty]]$circDensity
	allDensity	<- perTypeData$All$circDensity
	vertDensity	<- perTypeData$Ver$'3'$circDensity

	# Number of genes for each density plot
	numGenes	<- length(lHGTDensity$data)

	# Produce plot
	position_plot	<- circularDensityPlot(
		dataDensityA = lHGTDensity,
		dataDensityB = vertDensity,
		bgDensity = allDensity,
		shrink = 0.7,
		tcl.offset = 0.8,
		titleName = paste0("lHGT Enrichment at Penalty = ", penalty, "\nGenes = ", numGenes),
		densBLineColor = dataTypeCols$Ver)
	replayPlot(position_plot)
}))

outputFileName	<- file.path(circDensityFig_path, "lHGTbyPenalty.pdf")
invisible(dev.copy2pdf(file = outputFileName))
invisible(dev.off())

message("\rPlotting lHGT gene density across penalties... done")

# ------------------------------------------------------------------------------------- #
# Compare lHGT and sHGT density at a given penalty
message("Plotting lHGT vs sHGT gene density at single penatly...", appendLF = FALSE)

penalty	<- "4"

quartz(width = 12, height = 6.5)
par(mfrow = c(1, 2))
par(mar = c(0, 0, 0, 0))
invisible(lapply(list("lHGT", "sHGT"), function(dataType) {

	# Density values needed
	HGTDensity	<- perTypeData[[dataType]][[penalty]]$circDensity
	allDensity	<- perTypeData$All$circDensity
	vertDensity	<- perTypeData$Ver$'3'$circDensity

	# Number of genes for each density plot
	numGenes	<- length(HGTDensity$data)

	# Produce plot
	position_plot	<- circularDensityPlot(
		dataDensityA = HGTDensity,
		dataDensityB = vertDensity,
		bgDensity = allDensity,
		shrink = 0.7,
		tcl.offset = 0.8,
		titleName = paste0(dataType, " Enrichment at Penalty = ",penalty, "\nGenes = ", numGenes),
		densBLineColor = dataTypeCols$Ver)
	replayPlot(position_plot)
}))

outputFileName	<- file.path(circDensityFig_path, "lHGTvsHGT.pdf")
invisible(dev.copy2pdf(file = outputFileName))
invisible(dev.off())

message("\rPlotting lHGT vs sHGT gene density at single penatly... done")

# ------------------------------------------------------------------------------------- #
# Compare group vs subgroup lHGT density at a given penalty
# message("Plotting Old vs Recent lHGT gene density...", appendLF = FALSE)

# penalty		<- "4"
# # Define the bandwidth for circular density estimation
# bandwidth 	<- 3000

# quartz(width = 21, height = 8)
# par(mfrow = c(1, 3))
# par(mar = c(0, 0, 0, 0))
# invisible(lapply(list("lHGT", "Group", "Subgroup"), function(type) {

# 	# Determine the data set required (densities for Group / Subgroup have to be calculated)
# 	dataShortcut	<- perTypeData$lHGT[[penalty]]$allPosData
# 	if (identical(type, "Group")) {
# 		typeCircStart	<- dataShortcut$CircStart[which(dataShortcut$Subgroup == FALSE)]
# 	} else if (identical(type, "Subgroup")) {
# 		typeCircStart	<- dataShortcut$CircStart[which(dataShortcut$Subgroup == TRUE)]
# 	} else {
# 		typeCircStart	<- dataShortcut$CircStart
# 	}

# 	# Calculate circular density (on the start position)
# 	HGTDensity	<- density.circular(typeCircStart, kernel = "vonmises", bw = bandwidth)
# 	allDensity	<- perTypeData$All$circDensity
# 	vertDensity	<- perTypeData$Ver$'3'$circDensity

# 	# Number of genes for each density plot
# 	numGenes	<- length(HGTDensity$data)

# 	# Produce plot
# 	position_plot	<- circularDensityPlot(
# 		dataDensityA = HGTDensity,
# 		dataDensityB = vertDensity,
# 		bgDensity = allDensity,
# 		shrink = 0.7,
# 		tcl.offset = 0.8,
# 		titleCex = 1.8,
# 		titleName = paste0(type, " Enrichment at Penalty = ", penalty, "\nGenes = ", numGenes),
# 		densBLineColor = dataTypeCols$Ver)
# 	replayPlot(position_plot)
# }))

# outputFileName	<- file.path(circDensityFig_path, "OldvRecent_lHGT.pdf")
# invisible(dev.copy2pdf(file = outputFileName))
# invisible(dev.off())

# message("\rPlotting Old vs Recent lHGT gene density... done")

# message(paste0("\nAll circular density plots written to: ", circDensityFig_path, "\n"))

# ------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------- #
# # Make a plot to show how the subcompartments (subdivisions) overlap with the HGT density around genome
# message("Plotting lHGT gene density over subcompartment divisions...", appendLF = FALSE)

# subDivisionKey_df	<- subDivisionKey_data$subDivisionKey_df
# subDivison_cols		<- subDivisionKey_data$subDivison_cols

# # Prepare a bare circular plot with the subcompartment zones
# subDivisionKeyBare_plot		<- ggplot(data = subDivisionKey_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = col, col = NA)) +
# 	coord_polar("x") +
# 	scale_x_continuous(labels = NULL, breaks = c(0, pi), limits = c(0, 2 * pi)) +
# 	scale_y_continuous(limits = c(0, 1)) +
# 	geom_rect(col = NA, size = 0) +
# 	scale_fill_manual(values = alpha(subDivison_cols, 0.3), guide = FALSE) +
# 	darkTheme +
# 	theme(
# 		panel.grid = element_blank(),
# 		panel.grid.major = element_blank(),
# 		panel.grid.minor = element_blank(),
# 		axis.ticks = element_blank(),
# 		axis.text.y = element_blank(),
# 		axis.text.x = element_text(size = 18),
# 		panel.border = element_blank(),
# 		legend.key.size = unit(1, "cm"),
# 		legend.title = element_blank()
# 	)

# # Get densities, use penalty 4 for lHGT genes
# penalty		<- "4"

# lHGTDensity	<- perTypeData$lHGT[[penalty]]$circDensity
# vertDensity	<- perTypeData$Ver$'3'$circDensity
# allDensity	<- perTypeData$All$circDensity

# # Number of genes for each density plot
# numGenes	<- length(lHGTDensity$data)

# # Plotting this is tricky: straight-to-pdf gives split output. Plot to quartz device and read off device to PDF file
# quartz(width = 9, height = 9)
# plot.new()
# print(subDivisionKeyBare_plot)

# par(mar = c(0,0,0,0))
# par(new = TRUE)

# circularDensityPlot(dataDensityA = lHGTDensity, dataDensityB = vertDensity, bgDensity = allDensity, shrink = 0.9, tcl.offset = 0.8, titleName = paste0("Vertical Enrichment at Penalty = ", penalty, "\nGenes = ", numGenes), bg = NA)

# outputFileName	<- file.path(circDensityFig_path, "DensityOverSubcompartment.pdf")
# invisible(dev.copy2pdf(file = outputFileName))
# invisible(dev.off())

# message("\rPlotting lHGT gene density over subcompartment divisions... done")

# ------------------------------------------------------------------------------------- #
# # Warning message because file needs to be edited

# if (file.exists(file.path(figureOutput_path, "HGT_density_bySubdivision.ai"))) {
# 	message(paste0("\nThe HGT density plotted over subcompartments has been edited to produce:\n\t", file.path(figureOutput_path, "HGT_density_bySubdivision.ai")))
# } else {
# 	message(paste0("\nWARNING: The HGT density over subcompartment plot \"DensityOverSubcompartment.pdf\" should be adjusted / edited and placed in:\n\t", figureOutput_path))
# }



# ------------------------------------------------------------------------------------- #

binomial_list		<- unique(perTypeData$All$allPosData$binomial)

if (is.null(outlierTaxid)) {
	uinSpecies			<- 0.8
	titleCexSpecies		<- 0.7
	bandwidthSpecies	<- 1500
	shrinkSpecies		<- 1

	quartz(width = 12, height = 12)
	par(mfrow = c(5, 5))
	par(mar = c(0, 0, 0, 0))
	invisible(lapply(binomial_list, function(species) {

		# HGT and background circular data
		HGTSpec_circ	<- subset(perTypeData$lHGT$'4'$allPosData, binomial == species, select = CircStart, drop = TRUE)
		bgSpec_circ		<- subset(perTypeData$All$allPosData, binomial == species, select = CircStart, drop = TRUE)

		# HGT and background density calculation. NB bandwidth is lower than for the combined plots!
		HGTSpec_dens	<- density.circular(HGTSpec_circ, kernel = "vonmises", bw = bandwidthSpecies)
		bgSpec_dens		<- density.circular(bgSpec_circ, kernel = "vonmises", bw = bandwidthSpecies)

		# Number of genes for each density plot
		numGenes		<- length(HGTSpec_circ)

		# Produce plot
		position_plot	<- circularDensityPlot(dataDensityA = HGTSpec_dens, bgDensity = bgSpec_dens, shrink = shrinkSpecies, tcl.offset = 0.8, titleCex = titleCexSpecies, uin = uinSpecies, titleName = paste0(species, "\nGenes = ", numGenes))
		replayPlot(position_plot)
	}))
	outputFileName	<- file.path(circDensityFig_path, "bySpecies_HGT_Enrichment.pdf")
	invisible(dev.copy2pdf(file = outputFileName))
	invisible(dev.off())


	quartz(width = 12, height = 12)
	par(mfrow = c(5, 5))
	par(mar = c(0, 0, 0, 0))
	invisible(lapply(binomial_list, function(species) {

		# HGT and background circular data
		VerSpec_circ	<- subset(perTypeData$Ver$'3'$allPosData, binomial == species, select = CircStart, drop = TRUE)
		bgSpec_circ		<- subset(perTypeData$All$allPosData, binomial == species, select = CircStart, drop = TRUE)

		# HGT and background density calculation. NB bandwidth is lower than for the combined plots!
		VerSpec_dens	<- density.circular(VerSpec_circ, kernel = "vonmises", bw = bandwidthSpecies)
		bgSpec_dens		<- density.circular(bgSpec_circ, kernel = "vonmises", bw = bandwidthSpecies)

		# Number of genes for each density plot
		numGenes		<- length(VerSpec_circ)

		# Produce plot
		position_plot	<- circularDensityPlot(
			dataDensityA = VerSpec_dens,
			bgDensity = bgSpec_dens,
			enrichUpColor = wes_palette("BottleRocket2")[1],
			enrichDownColor = wes_palette("Zissou1")[1],
			shrink = shrinkSpecies,
			tcl.offset = 0.8,
			titleCex = titleCexSpecies,
			uin = uinSpecies,
			titleName = paste0(species, "\nGenes = ", numGenes))
		replayPlot(position_plot)
	}))
	outputFileName	<- file.path(circDensityFig_path, "bySpecies_Ver_Enrichment.pdf")
	invisible(dev.copy2pdf(file = outputFileName))
	invisible(dev.off())
}



# # Just G kaustophilus vertical density
# quartz(width = 14, height = 14, canvas = "white", bg = "white")
# par(mfrow = c(1, 1))
# par(mar = c(0, 0, 0, 0))

# invisible(lapply(binomial_list[1], function(species) {

# 	# HGT and background circular data
# 	VerSpec_circ	<- subset(perTypeData$Ver$'3'$allPosData, binomial == species, select = CircStart, drop = TRUE)
# 	bgSpec_circ		<- subset(perTypeData$All$allPosData, binomial == species, select = CircStart, drop = TRUE)

# 	# HGT and background density calculation. NB bandwidth is lower than for the combined plots!
# 	VerSpec_dens	<- density.circular(VerSpec_circ, kernel = "vonmises", bw = 500)
# 	bgSpec_dens		<- density.circular(bgSpec_circ, kernel = "vonmises", bw = 500)

# 	# Number of genes for each density plot
# 	numGenes		<- length(VerSpec_circ)

# 	# Produce plot
# 	position_plot	<- circularDensityPlot(
# 		dataDensityA = VerSpec_dens,
# 		bgDensity = bgSpec_dens,
# 		bg = NULL,
# 		axisCol = wes_palette("Moonrise1")[4],
# 		enrichUpColor = wes_palette("Darjeeling1")[3],
# 		enrichDownColor = wes_palette("FantasticFox1")[3],
# 		shrink = 0.7,
# 		tcl.offset = 0.9,
# 		titleCex = 1.5,
# 		uin = 3.5,
# 		axis.at = seq(0, 359, by = 30),
# 		axis.labels = seq(0, 359, by = 30),
# 		titleName = paste0(species, "\nGenes = ", numGenes))
# 	replayPlot(position_plot)
# }))

# outputFileName	<- file.path(circDensityFig_path, "gKaustophlis_verticalGenes")
# invisible(quartz.save(file = paste0(outputFileName, ".pdf"), type = "pdf", dpi = 100))
# invisible(dev.off())




