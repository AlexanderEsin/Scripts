#!/usr/bin/env Rscript

# Load master variables and HGT position functions
invisible(sapply(HGTPos.all, source, .GlobalEnv))

require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load("dplyr", "wesanderson", "ggplot2")

# ------------------------------------------------------------------------------------- #
# Read in data

message("\nReading in data...", appendLF = FALSE)

# General position data
perTypeData			<- readRDS(file.path(positionData_path, "AG_perTypeData.rds"))

# Subgroup data
subgroupData		<- readRDS(file.path(positionData_path, "AG_subgroupData.rds"))

# Subdivision key
subDivisionKey_file	<- file.path(positionData_path, "subDivisionKeyData.rds")
if(!file.exists(subDivisionKey_file)) {
	stop("The subdivision data file \"", subDivisionKey_file, "\" is required.\nRun the subDivisionPrep.R script!")
} else {
	subDivisionKey_data	<- readRDS(file.path(positionData_path, "subDivisionKeyData.rds"))
}

message("\rReading in data... done\n")

# ------------------------------------------------------------------------------------- #
# Plotting

# Output path for figures
SBiasFigPath	<- file.path(figureOutput_path, "Strand_bias")

# Quartz plotting options common to this script
quartz.options(canvas = "#333233", bg = "#333233", type = "png", dpi = 300)

# ------------------------------------------------------------------------------------- #
# Convert the compartment key to linear

xValCols		<- subset(subDivisionKey_data$subDivisionKey_df, select = c(xmin, xmax))
xValLineCols	<- 100 * (xValCols / (2 * pi))
subDivisionKeyLinear_df	<- cbind(xValLineCols, subset(subDivisionKey_data$subDivisionKey_df, select = -c(xmin, xmax)))

## For the strandBias plots, the bars will be drawn across most of the scale
subDivisionKeyLinear_df$ymin	<- -1
subDivisionKeyLinear_df$ymax	<- 1

# Rename the subcompartments
levels(subDivisionKeyLinear_df$col)	<- c("Ori", "Near Ori", "Far Ori", "Far Ter", "Near Ter", "Ter")

# Also get the subdivision colors
subDivison_cols	<- subDivisionKey_data$subDivison_cols

# ------------------------------------------------------------------------------------- #

binomial_list			<- unique(perTypeData$All$allPosData$binomial)

# Calculate strand biad per type of gene (including Old and Recent lHGTS) and across every genome (binomial)
byTypeStrandBias_list	<- lapply(dataTypes_withAge, processStrandBias, data = perTypeData)
names(byTypeStrandBias_list)	<- dataTypes_withAge

# Extract the roll mean values (for outlier ID), and then adjust the All dataframe to remove the extra column
outlierID_rollmean			<- byTypeStrandBias_list$All
byTypeStrandBias_list$All	<- byTypeStrandBias_list$All[, !(names(byTypeStrandBias_list$All) == "rollMean")]

# Subset outliers
outlierSpecies		<- subset(outlierID_rollmean, BinIndex == 50 & (rollMean > 0.35 | rollMean < -0.35), select = Species, drop = TRUE)
outlierLineTypes	<- c("dashed", "twodash")
if (length(outlierSpecies) != length(outlierLineTypes)) {
	stop(paste0("Change number of linetypes available to: ", length(outlierSpecies)))
}

# Define a new column - in which we only have the outlier Species names (to use as linetype aesthetic below)
byTypeStrandBias_list	<- lapply(byTypeStrandBias_list, function(dataType) {
	dataType$lineType	<- NA
	dataType$lineType[which(dataType$Species %in% outlierSpecies)]	<- dataType$Species[which(dataType$Species %in% outlierSpecies)]
	return(dataType)
})

# ------------------------------------------------------------------------------------- #
# Prepare plot data and color palettes

SB_allGenes			<- byTypeStrandBias_list$All
SB_allGenes$Type	<- factor(SB_allGenes$Type, levels = dataTypes_withAge)

SB_AllvsVer			<- bind_rows(list(byTypeStrandBias_list$All, byTypeStrandBias_list$Ver))
SB_AllvsVer$Type	<- factor(SB_AllvsVer$Type, levels = dataTypes_withAge)

SB_VervsHGT			<- bind_rows(list(byTypeStrandBias_list$Ver, byTypeStrandBias_list$lHGT))
SB_VervsHGT$Type	<- factor(SB_VervsHGT$Type, levels = dataTypes_withAge)

SB_OldvsNew			<- bind_rows(list(byTypeStrandBias_list$Old, byTypeStrandBias_list$Recent))
SB_OldvsNew$Type	<- factor(SB_OldvsNew$Type, levels = dataTypes_withAge)


allGenes_col	<- "#D9D9D9"
allVsVer_cols	<- c("#D9D9D9", wes_palette("Darjeeling1")[2])
verVsHGT_cols	<- wes_palette("Darjeeling1")[1:2]
oldVsNew_cols	<- wes_palette("Darjeeling1")[4:5]

# ------------------------------------------------------------------------------------- #
### Produce plots

## Plot strand bias across all genes: highlight a couple of outliers. Need to be cross-referenced with Mauve output
SB_AllGenes_plot	<- ggplot(SB_allGenes, aes(x = BinIndex, y = StrandBias, color = Type)) +
	scale_x_continuous(name = "Bin along genome", limits = c(0, 100), breaks = seq(0, 100, 10), minor_breaks = NULL) +
	scale_y_continuous(name = "Enrichment for genes on same strand as DnaA") +
	geom_hline(yintercept = 0, col = alpha("#D9D9D9", 0.75), size = 1) +
	geom_rect(data =  subDivisionKeyLinear_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = alpha(unlist(lapply(subDivison_cols, rep, times = 2)), 0.2), inherit.aes = FALSE) +
	# Most species smooth lines
	geom_smooth(
		data = subset(SB_allGenes, !(Species %in% outlierSpecies)),
		aes(group = Species),
		method = "loess",
		span = 0.25,
		col = alpha(allGenes_col, 0.5),
		size = 0.2,
		se = FALSE) +
	# Outlier species smooth lines
	geom_smooth(
		data = subset(SB_allGenes, (Species %in% outlierSpecies)),
		aes(linetype = lineType),
		method = "loess",
		span = 0.25,
		col = alpha(allGenes_col, 0.5),
		size = 0.5,
		se = FALSE) +
	# Overall smooth line (cross-species)
	geom_smooth(aes(color = Type, fill = Type), method = "loess", size = 1.5, se = TRUE, span = 0.25) +
	scale_color_manual(values = allGenes_col, guide = guide_legend(title = "Gene Type")) +
	scale_fill_manual(values = alpha(allGenes_col, 0.7), guide = FALSE) +
	scale_linetype_manual(values = outlierLineTypes, guide = guide_legend(title = "Outliers")) +
	darkTheme +
	theme(
		legend.justification = c(1, 1),
		legend.position = c(1, 0.825),
		legend.background = element_rect(fill = "#333233", color = "#D9D9D9")
	)


## Compare All genes to Vertical genes
SB_AllvsVer_plot	<- ggplot(SB_AllvsVer, aes(x = BinIndex, y = StrandBias, color = Type)) +
	scale_x_continuous(name = "Bin along genome", limits = c(0, 100), breaks = seq(0, 100, 10), minor_breaks = NULL) +
	scale_y_continuous(name = "Enrichment for genes on same strand as DnaA") +
	geom_hline(yintercept = 0, col = alpha("#D9D9D9", 0.75), size = 1) +
	geom_rect(data =  subDivisionKeyLinear_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = alpha(unlist(lapply(subDivison_cols, rep, times = 2)), 0.2), inherit.aes = FALSE) +
	# Smooth lines for "All" genes: per species
	geom_smooth(
		data = subset(SB_AllvsVer, Type == "All"),
		aes(group = Species),
		method = "loess",
		span = 0.25,
		col = alpha(allVsVer_cols[1], 0.7),
		size = 0.2,
		se = FALSE) +
	# Smooth lines for "Ver" genes: per species
	geom_smooth(
		data = subset(SB_AllvsVer, Type == "Ver"),
		aes(group = Species),
		method = "loess",
		span = 0.25,
		col = alpha(allVsVer_cols[2], 0.7),
		size = 0.2,
		se = FALSE) +
	# Smooth lines for all Types of genes: cross species
	geom_smooth(aes(color = Type, fill = Type), method = "loess", size = 1.5, se = TRUE, span = 0.25) +
	scale_color_manual(values = allVsVer_cols, guide = guide_legend(title = "Gene Type")) +
	scale_fill_manual(values = alpha(allVsVer_cols, 0.7), guide = FALSE) +
	darkTheme +
	theme(
		legend.justification = c(1, 1),
		legend.position = c(1, 0.825),
		legend.background = element_rect(fill = "#333233", color = "#D9D9D9")
	)


## Compare Vertical genes to lHGT genes
SB_VervsHGT_plot	<- ggplot(SB_VervsHGT, aes(x = BinIndex, y = StrandBias, color = Type)) +
	scale_x_continuous(name = "Bin along genome", limits = c(0, 100), breaks = seq(0, 100, 10), minor_breaks = NULL) +
	scale_y_continuous(name = "Enrichment for genes on same strand as DnaA") +
	geom_hline(yintercept = 0, col = alpha("#D9D9D9", 0.75), size = 1) + 
	geom_rect(data =  subDivisionKeyLinear_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = alpha(unlist(lapply(subDivison_cols, rep, times = 2)), 0.2), inherit.aes = FALSE) +
	# Smooth lines for "lHGT" genes: per species
	geom_smooth(
		data = subset(SB_VervsHGT, Type == "lHGT"),
		aes(group = Species),
		method = "loess",
		span = 0.25,
		color = alpha(verVsHGT_cols[1], 0.8),
		size = 0.2,
		se = FALSE) +
	# Smooth lines for "Ver" genes: per species
	geom_smooth(
		data = subset(SB_VervsHGT, Type == "Ver"),
		aes(group = Species),
		method = "loess",
		span = 0.25,
		color = alpha(verVsHGT_cols[2], 0.8),
		size = 0.2,
		se = FALSE) +
	# Smooth lines for all Types of genes: cross species
	geom_smooth(aes(color = Type, fill = Type), method = "loess", size = 1.5, se = TRUE, span = 0.25) +
	scale_color_manual(values = verVsHGT_cols, guide = guide_legend(title = "Gene Type")) +
	scale_fill_manual(values = alpha(verVsHGT_cols, 0.7), guide = FALSE) + 
	darkTheme +
	theme(
		legend.justification = c(1, 1),
		legend.position = c(0.9, 0.9),
		legend.background = element_rect(fill = "#333233", color = "#D9D9D9")
	)


## Compare Old lHGTs to Recent lHGTS. NB: Recent lHGTS are sparse - hence large SE and inclusion of histogram 
SB_OldvsNew_plot	<- ggplot(SB_OldvsNew, aes(x = BinIndex, y = StrandBias, color = Type)) +
	scale_x_continuous(name = "Bin along genome", limits = c(0, 100), breaks = seq(0, 100, 10), minor_breaks = NULL) +
	scale_y_continuous(name = "Enrichment for genes on same strand as DnaA", limits = c(-1, 1.5)) +
	geom_rect(data =  subDivisionKeyLinear_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = alpha(unlist(lapply(subDivison_cols, rep, times = 2)), 0.2), inherit.aes = FALSE) +
	# Histogram for Old genes: show per-Bin population (relative)
	geom_histogram(
		data = subset(SB_OldvsNew, Type == "Old" & !(is.na(StrandBias))),
		aes(x = BinIndex,y = (..ncount.. / 4)),
		binwidth = 1,
		fill = alpha(oldVsNew_cols[1], 0.4),
		color = alpha("#333233", 0.8),
		size = 0.2,
		inherit.aes = FALSE) +
	# Histogram for Recent genes: show per-Bin population (relative)
	geom_histogram(
		data = subset(SB_OldvsNew, Type == "Recent" & !(is.na(StrandBias))),
		aes(x = BinIndex,y = -(..ncount.. / 4)),
		binwidth = 1,
		fill = alpha(oldVsNew_cols[2], 0.4),
		color = alpha("#333233", 0.8),
		size = 0.2,
		inherit.aes = FALSE) +
	# Plot 0-line later to go over bottom of histograms
	geom_hline(yintercept = 0, col = alpha("#D9D9D9", 0.75), size = 1) +
	# Smooth lines for all Types of genes: cross species
	geom_smooth(aes(color = Type, fill = Type), method = "loess", size = 1.5, span = 0.25) +
	scale_color_manual(values = oldVsNew_cols, guide = guide_legend(title = "Gene Type")) +
	scale_fill_manual(values = alpha(oldVsNew_cols, 0.7), guide = FALSE) + 
	darkTheme +
	theme(
		legend.justification = c(1, 1),
		legend.position = c(0.9, 0.9),
		legend.background = element_rect(fill = "#333233", color = "#D9D9D9")
	)

# ------------------------------------------------------------------------------------- #
# Write plots

quartz(width = 16, height = 8, file = file.path(SBiasFigPath, "All_genes.png"))
print(SB_AllGenes_plot)
invisible(dev.off())

quartz(width = 16, height = 8, file = file.path(SBiasFigPath, "AllvsVertical.png"))
print(SB_AllvsVer_plot)
invisible(dev.off())

quartz(width = 16, height = 8, file = file.path(SBiasFigPath, "VerticalvslHGT.png"))
print(SB_VervsHGT_plot)
invisible(dev.off())

quartz(width = 16, height = 8, file = file.path(SBiasFigPath, "OldvsRecent_HGT.png"))
print(SB_OldvsNew_plot)
invisible(dev.off())









