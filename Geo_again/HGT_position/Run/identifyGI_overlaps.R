#!/usr/bin/env Rscript

# Load master variables and HGT position functions
invisible(sapply(HGTPos.all, source, .GlobalEnv))

require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load("dplyr", "ggplot2", "wesanderson")

# ------------------------------------------------------------------------------------- #`
# Read in data

message("\nReading in data...", appendLF = FALSE)

# General position data
perTypeData			<- readRDS(file.path(positionData_path, "AG_perTypeData.rds"))

# GI-specific data
byType_withGI_data	<- readRDS(file.path(positionData_path, "AG_perTypeGIData.rds"))
GI_positions_data	<- readRDS(file.path(positionData_path, "AG_GI_positions.rds"))

message("\rReading in data... done\n")

# ------------------------------------------------------------------------------------- #
# Output path for GI figures
GIanalysisFig_path	<- file.path(figureOutput_path, "genomicIslands")
if (!dir.exists(GIanalysisFig_path)) dir.create(GIanalysisFig_path)

# ------------------------------------------------------------------------------------- #
# Overall density distribution of genomic islands across genomes
GI_positions_parts	<- GI_positions_data$GI_positions_data

GI_crossGenomeDens_plot	<- ggplot(data = GI_positions_parts, aes(x = GI_location)) +
	scale_x_continuous(
		name = "Relative Genomic Position",
		limits = c(0, 1),
		breaks = seq(0, 1, by = 0.5),
		minor_breaks = seq(0, 1, by = ((1 / 360) * 30)),
		labels = c("Origin", "Terminus", "Origin"),
		sec.axis = dup_axis(
			name = NULL,
			breaks = seq(0, 1, by = ((1 / 360) * 30)),
			labels = c("Origin", seq(30, 150, by = 30), "Terminus", seq(210, 330, by = 30) ,"Origin")
		)
	) +
	geom_histogram(aes(y = ..ncount.. / 4, fill = binomial), bins = 1000) +
	scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
	stat_density(geom = "line", adjust = 1/4, color = wes_palette("Moonrise2")[4], show.legend = FALSE) +
	scale_fill_manual(values = alpha(wes_palette("Darjeeling1", n = numSpeciesGI, type = "continuous"), 0.7)) +
	theme_classic() +
	theme(
		panel.grid.major.y = element_line(size = 0.05),
		panel.grid.minor.x  = element_line(size = 0.4, linetype = "dashed", color = wes_palette("Darjeeling1")[2]),
		axis.text = element_text(size = 11),
		axis.title = element_text(size = 14),
		axis.title.y = element_blank(),
		axis.ticks.y = element_blank(),
		axis.ticks.length = unit(0.6, "lines"),
		plot.title = element_text(hjust = 0.5, size = 16))


quartz(width = 20, height = 10, canvas = "white", bg = "white")
print(GI_crossGenomeDens_plot)
invisible(quartz.save(file = file.path(GIanalysisFig_path, "GI_crossGenomePositions.pdf"), type = "pdf", dpi = 100))
invisible(dev.off())

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
		subset_data	<- byType_withGI_data$lHGT[[lHGT_penalty]]
	} else {
		subset_data	<- subset(byType_withGI_data$lHGT[[lHGT_penalty]], In_GI == getInGI)
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
	subset_data	<- subset(byType_withGI_data$lHGT[[penalty]], In_GI == FALSE)

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































