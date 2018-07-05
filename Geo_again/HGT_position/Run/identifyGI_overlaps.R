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

# Zone list
zoneBoundaryList	<- readRDS(file.path(positionData_path, "AG_zoneBoundaries.rds"))

message("\rReading in data... done\n")

# ------------------------------------------------------------------------------------- #
# Output path for GI figures
GIanalysisFig_path	<- file.path(figureOutput_path, "genomicIslands")
if (!dir.exists(GIanalysisFig_path)) dir.create(GIanalysisFig_path)

# ------------------------------------------------------------------------------------- #
# Overall density distribution of genomic islands across genomes
GI_positions_parts	<- GI_positions_data$GI_positions_data
speciesWithGI		<- unique(GI_positions_data$GI_positions_data$binomial)
numSpeciesGI		<- length(speciesWithGI)

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
quartz(width = 21, height = 10)
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
	circularDensityPlot(
		dataDensityA = subset_dens,
		bgDensity = allPosData_withGI_dens,
		shrink = 0.7,
		tcl.offset = 0.8,
		titleCex = 1.8,
		bg = "white",
		axisCol = axisCol,
		enrichUpColor = dataTypeCols$HGT,
		enrichDownColor = wes_palette("Darjeeling1")[1],
		densBLineColor = dataTypeCols$Ver,
		titleName = paste0(subset, "\nGenes = ", nrow(subset_data)))
}))

outputFileName	<- file.path(GIanalysisFig_path, "GIvsNon-GI_lHGT_density.pdf")
invisible(dev.copy2pdf(file = outputFileName))
invisible(dev.off())


# ------------------------------------------------------------------------------------- #
# Compare the Outside GI lHGT gene distribution across multiple penalties


# Set up plot
quartz(width = 14, height = 14)
par(mfrow = c(2, 2))
par(mar = c(0, 0, 0, 0))

invisible(lapply(penalty_list, function(penalty) {

	# Subset data by penalty
	subset_data	<- subset(byType_withGI_data$lHGT[[penalty]], In_GI == FALSE)

	# Calculate circular density
	subset_dens	<- density.circular(subset_data$CircStart, kernel = "vonmises", bw = bandwith)

	# Produce circular plot
	circularDensityPlot(
		dataDensityA = subset_dens,
		bgDensity = allPosData_withGI_dens,
		shrink = 0.7,
		tcl.offset = 0.8,
		titleCex = 1.8,
		bg = "white",
		axisCol = axisCol,
		enrichUpColor = dataTypeCols$HGT,
		enrichDownColor = wes_palette("Darjeeling1")[1],
		densBLineColor = dataTypeCols$Ver,
		titleName = paste0("Penalty = ", penalty, "\nGenes = ", nrow(subset_data)))
}))

outputFileName	<- file.path(GIanalysisFig_path, "crossPenaltyNon-GI_density.pdf")
invisible(dev.copy2pdf(file = outputFileName))
invisible(dev.off())

# ------------------------------------------------------------------------------------- #

# Total HGT genes in the 17 GI-ready genomes
allHGT_genes	<- byType_withGI_data$lHGT$`4`

allVer_genes	<- byType_withGI_data$Ver$`3`

# Number of all genes within
allGI_genes		<- byType_withGI_data$HGT %>% subset(In_GI == TRUE)

# Overlap between GIs within HGT zones and GIs outside of GI zones
vennTemp	<- venn.diagram(
	list(HGT_genes = allHGT_genes$protID, GI_genes = allGI_genes$protID),
	filename = NULL,
	fill = c(dataTypeCols$HGT, "darkblue"),
	alpha = c(0.5, 0.5))

quartz(h = 8, w = 8)
grid.draw(vennTemp)
outputFileName	<- file.path(GIanalysisFig_path, "HGT_vs_GI_vennOverlap.pdf")
invisible(dev.copy2pdf(file = outputFileName))
invisible(dev.off())


# ------------------------------------------------ #

numGIs_byZone	<- zoneBoundaryList$fullRange %>%
	group_by(zoneMin) %>%
	mutate(zoneArea = zoneMax - zoneMin) %>%
	mutate(numGIs = nrow(allGI_genes[findInterval(allGI_genes$relGeneStart, c(zoneMin, zoneMax)) == 1L,]))

numGIs_byType	<- numGIs_byZone %>%
	group_by(zoneType) %>%
	summarise(GIsByType = sum(numGIs), areaByType = sum(zoneArea)) %>%
	mutate(GI_density = round(GIsByType / areaByType, 1))

byTypeBoundary_split	<- split(numGIs_byZone, numGIs_byZone$zoneType)


# ------------------------------------------------ #

withinSeqOfRanges	<- function(value, leftEnd, rightEnd) {
	if (length(leftEnd) != length(rightEnd)) {
		stop("The edges of the ranges have to be the same length")
	}

	rightmost.closed = FALSE
	left.open = FALSE

	valueWithinRange	<- lapply(1:length(leftEnd), function(rangeIndex) {
		leftEdge	<- leftEnd[rangeIndex]
		rightEdge	<- rightEnd[rangeIndex]

		if (rightEdge == 1) {
			rightmost.closed = TRUE
		}

		withinRange	<- which(findInterval(value, c(leftEdge, rightEdge), rightmost.closed = rightmost.closed) == 1)
		return(withinRange)
	})

	return(unlist(valueWithinRange))
}

allGenesByZone <- byType_withGI_data$All %>% mutate(zone = case_when(
	lapply(relGeneStart, withinSeqOfRanges, leftEnd = byTypeBoundary_split$HGT$zoneMin, rightEnd = byTypeBoundary_split$HGT$zoneMax) == 1 ~ "HGT",
	lapply(relGeneStart, withinSeqOfRanges, leftEnd = byTypeBoundary_split$Other$zoneMin, rightEnd = byTypeBoundary_split$Other$zoneMax) == 1 ~ "Other",
	lapply(relGeneStart, withinSeqOfRanges, leftEnd = byTypeBoundary_split$Ver$zoneMin, rightEnd = byTypeBoundary_split$Ver$zoneMax) == 1 ~ "Ver"
	)
)

HGTzone_allGenes	<- allGenesByZone %>% subset(zone == "HGT")
Otherzone_allGenes	<- allGenesByZone %>% subset(zone == "Other")


vennTemp	<- venn.diagram(
	list(HGT_zones = HGTzone_allGenes$protID, GI_genes = allGI_genes$protID, Other_zones = Otherzone_allGenes$protID),
	filename = NULL,
	fill = c(dataTypeCols$HGT, "blue", dataTypeCols$Other),
	alpha = c(0.5, 0.5, 0.5))

quartz(h = 8, w = 8)
grid.draw(vennTemp)
outputFileName	<- file.path(GIanalysisFig_path, "HGTzone_OtherZone_vs_GI_vennOverlap.pdf")
invisible(dev.copy2pdf(file = outputFileName))
invisible(dev.off())


# ------------------------------------------------------------------------------------- #
# Ran a test with the downloaded GI data for GKau (with locus IDs from Island Viewer)
# They use an earlier assembly revision, so 4 proteins we have they are lacking. Thus, 93/97 proteins overlap

gkauGIViewer_data	<- read.table(file.path(giProcess_path, "testGKau.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# lHGTinGI_crossCheck	<- testlHGT_data[which(testlHGT_data$locusTag %in% gkauGIViewer_data$Locus),]
# lHGTinGI_myMeth[which(!lHGTinGI_myMeth$locusTag %in% lHGTinGI_crossCheck$locusTag),]































