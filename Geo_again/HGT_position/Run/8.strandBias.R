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

# Zone list
zoneBoundaryList	<- readRDS(file.path(positionData_path, "AG_zoneBoundaries.rds"))

message("\rReading in data... done\n")

# ------------------------------------------------------------------------------------- #
# Plotting

# Output path for figures
SBiasFigPath	<- file.path(figureOutput_path, "Strand_bias")
if (!dir.exists(SBiasFigPath)) dir.create(SBiasFigPath)

# Quartz plotting options common to this script
quartz.options(canvas = "white", bg = "white")

# ------------------------------------------------------------------------------------- #

binomial_list	<- unique(perTypeData$All$allPosData$binomial)
binNumber		<- 200

# Calculate strand biad per type of gene (including Old and Recent lHGTS) and across every genome (binomial)
byTypeStrandBias_list			<- lapply(dataTypes_withAge, processStrandBias, data = perTypeData, binNumber = binNumber)
names(byTypeStrandBias_list)	<- dataTypes_withAge

# Extract the roll mean values (for outlier ID), and then adjust the All dataframe to remove the extra column
outlierID_rollmean			<- byTypeStrandBias_list$All
byTypeStrandBias_list$All	<- byTypeStrandBias_list$All[, !(names(byTypeStrandBias_list$All) == "rollMean")]

# Subset outliers
outlierSpecies		<- subset(outlierID_rollmean, BinIndex == 50 & (rollMean > 0.35 | rollMean < -0.35), select = Species, drop = TRUE)

#--There are no outlier species--#

# ------------------------------------------------------------------------------------- #
# Prepare plot data and color palettes

SB_allGenes			<- byTypeStrandBias_list$All
SB_allGenes_ext		<- rbind(SB_allGenes, cbind(SB_allGenes[which(SB_allGenes$BinIndex > 0 & SB_allGenes$BinIndex <= 0.25), 1, drop = FALSE] + 1, SB_allGenes[which(SB_allGenes$BinIndex > 0 & SB_allGenes$BinIndex <= 0.25), -1]))
SB_allGenes_ext$Type	<- factor(SB_allGenes_ext$Type, levels = dataTypes_withAge)

SB_AllvsVer			<- bind_rows(list(byTypeStrandBias_list$All, byTypeStrandBias_list$Ver))
SB_AllvsVer_ext		<- rbind(SB_AllvsVer, cbind(SB_AllvsVer[which(SB_AllvsVer$BinIndex > 0 & SB_AllvsVer$BinIndex <= 0.25), 1, drop = FALSE] + 1, SB_AllvsVer[which(SB_AllvsVer$BinIndex > 0 & SB_AllvsVer$BinIndex <= 0.25), -1]))
SB_allGenes_ext$Type	<- factor(SB_allGenes_ext$Type, levels = dataTypes_withAge)

SB_VervsHGT			<- bind_rows(list(byTypeStrandBias_list$Ver, byTypeStrandBias_list$lHGT))
SB_VervsHGT_ext		<- rbind(SB_VervsHGT, cbind(SB_VervsHGT[which(SB_VervsHGT$BinIndex > 0 & SB_VervsHGT$BinIndex <= 0.25), 1, drop = FALSE] + 1, SB_VervsHGT[which(SB_VervsHGT$BinIndex > 0 & SB_VervsHGT$BinIndex <= 0.25), -1]))
SB_allGenes_ext$Type	<- factor(SB_allGenes_ext$Type, levels = dataTypes_withAge)

SB_OldvsNew			<- bind_rows(list(byTypeStrandBias_list$Old, byTypeStrandBias_list$Recent))
SB_OldvsNew_ext		<- rbind(SB_OldvsNew, cbind(SB_OldvsNew[which(SB_OldvsNew$BinIndex > 0 & SB_OldvsNew$BinIndex <= 0.25), 1, drop = FALSE] + 1, SB_OldvsNew[which(SB_OldvsNew$BinIndex > 0 & SB_OldvsNew$BinIndex <= 0.25), -1]))
SB_allGenes_ext$Type	<- factor(SB_allGenes_ext$Type, levels = dataTypes_withAge)

# ------------------------------------------------------------------------------------- #
# Plot strand bias across all genes
SB_AllGenes_plot	<- ggplot(data = SB_allGenes_ext, aes(x = BinIndex, y = StrandBias, color = Type)) +
	scale_x_continuous(
		expand = expand_scale(mult = c(0.025, 0.025)),
		name = "Normalized genome position",
		limits = c(0, 1.25),
		breaks = seq(0, 1.25, by = 0.5),
		labels = c("Origin", "Terminus", "Origin"),
		minor_breaks = seq(0, 1.25, by = ((1 / 360) * 30)),
		sec.axis = dup_axis(
			name = NULL,
			breaks = rowMeans(zoneBoundaryList$expandedRange[c("zoneMin", "zoneMax")]),
			labels = zoneBoundaryList$expandedRange$zoneName)) +
	scale_y_continuous(name = "Enrichment for genes on same strand as DnaA") +
	geom_hline(yintercept = 0, col = alpha(axisCol, 0.75), size = 1) +
	geom_vline(xintercept = zoneBoundaryList$expandedRange$boundary, color = zoneBoundaryList$boundCol, linetype = "dashed") +	
	geom_rect(data = zoneBoundaryList$expandedRange, aes(xmin = zoneMin, xmax = zoneMax, ymin = -Inf, ymax = Inf), fill = zoneBoundaryList$expandedRange$zoneCol_alpha, inherit.aes = FALSE) +
	geom_smooth(
		data = SB_allGenes_ext,
		aes(group = Species),
		method = "loess",
		span = 0.1,
		col = alpha(dataTypeCols$All, 0.5),
		size = 0.2,
		se = FALSE) +
	geom_smooth(aes(color = Type, fill = Type), method = "loess", size = 1.5, se = TRUE, span = 0.1) +
	scale_color_manual(values = dataTypeCols$All, guide = guide_legend(title = "Gene Type")) +
	scale_fill_manual(values = alpha(dataTypeCols$All, 0.01), guide = FALSE) +
	guides(color = guide_legend(override.aes = list(fill = NA))) +
	lightTheme +
	theme(
		legend.justification = c(1, 1),
		legend.position = c(0.7, 0.825),
		legend.background = element_rect(fill = "white", color = axisCol)
	)


# ------------------------------------------------ #
# Compare All genes to Vertical genes
SB_AllvsVer_plot	<- ggplot(SB_AllvsVer_ext, aes(x = BinIndex, y = StrandBias, color = Type)) +
	scale_x_continuous(
		expand = expand_scale(mult = c(0.025, 0.025)),
		name = "Normalized genome position",
		limits = c(0, 1.25),
		breaks = seq(0, 1.25, by = 0.5),
		labels = c("Origin", "Terminus", "Origin"),
		minor_breaks = seq(0, 1.25, by = ((1 / 360) * 30)),
		sec.axis = dup_axis(
			name = NULL,
			breaks = rowMeans(zoneBoundaryList$expandedRange[c("zoneMin", "zoneMax")]),
			labels = zoneBoundaryList$expandedRange$zoneName)) +
	scale_y_continuous(name = "Enrichment for genes on same strand as DnaA") +
	geom_hline(yintercept = 0, col = alpha(axisCol, 0.75), size = 1) +
	geom_vline(xintercept = zoneBoundaryList$expandedRange$boundary, color = zoneBoundaryList$boundCol, linetype = "dashed") +	
	geom_rect(data = zoneBoundaryList$expandedRange, aes(xmin = zoneMin, xmax = zoneMax, ymin = -Inf, ymax = Inf), fill = zoneBoundaryList$expandedRange$zoneCol_alpha, inherit.aes = FALSE) +
	geom_smooth(
		data = subset(SB_AllvsVer_ext, Type == "All"),
		aes(group = Species),
		method = "loess",
		span = 0.1,
		col = alpha(dataTypeCols$All, 0.7),
		size = 0.2,
		se = FALSE) +
	# Smooth lines for "Ver" genes: per species
	geom_smooth(
		data = subset(SB_AllvsVer_ext, Type == "Ver"),
		aes(group = Species),
		method = "loess",
		span = 0.1,
		col = alpha(dataTypeCols$Ver, 0.7),
		size = 0.2,
		se = FALSE) +
	# Smooth lines for all Types of genes: cross species
	geom_smooth(aes(color = Type, fill = Type), method = "loess", size = 1.5, se = TRUE, span = 0.1) +
	scale_color_manual(values = c(dataTypeCols$All, dataTypeCols$Ver), guide = guide_legend(title = "Gene Type")) +
	scale_fill_manual(values = c(dataTypeCols$All, dataTypeCols$Ver), guide = FALSE) +
	guides(color = guide_legend(override.aes = list(fill = NA))) +
	lightTheme +
	theme(
		legend.justification = c(1, 1),
		legend.position = c(0.7, 0.825),
		legend.background = element_rect(fill = "white", color = axisCol)
	)


# ------------------------------------------------ #
# Compare Vertical genes to lHGT genes
SB_VervsHGT_plot	<- ggplot(SB_VervsHGT_ext, aes(x = BinIndex, y = StrandBias, color = Type)) +
	scale_x_continuous(
		expand = expand_scale(mult = c(0.025, 0.025)),
		name = "Normalized genome position",
		limits = c(0, 1.25),
		breaks = seq(0, 1.25, by = 0.5),
		labels = c("Origin", "Terminus", "Origin"),
		minor_breaks = seq(0, 1.25, by = ((1 / 360) * 30)),
		sec.axis = dup_axis(
			name = NULL,
			breaks = rowMeans(zoneBoundaryList$expandedRange[c("zoneMin", "zoneMax")]),
			labels = zoneBoundaryList$expandedRange$zoneName)) +
	scale_y_continuous(name = "Enrichment for genes on same strand as DnaA") +
	geom_hline(yintercept = 0, col = alpha(axisCol, 0.75), size = 1) +
	geom_vline(xintercept = zoneBoundaryList$expandedRange$boundary, color = zoneBoundaryList$boundCol, linetype = "dashed") +	
	geom_rect(data = zoneBoundaryList$expandedRange, aes(xmin = zoneMin, xmax = zoneMax, ymin = -Inf, ymax = Inf), fill = zoneBoundaryList$expandedRange$zoneCol_alpha, inherit.aes = FALSE) +
	geom_smooth(
		data = subset(SB_VervsHGT_ext, Type == "lHGT"),
		aes(group = Species),
		method = "loess",
		span = 0.1,
		color = alpha(dataTypeCols$HGT, 0.8),
		size = 0.2,
		se = FALSE) +
	# Smooth lines for "Ver" genes: per species
	geom_smooth(
		data = subset(SB_VervsHGT_ext, Type == "Ver"),
		aes(group = Species),
		method = "loess",
		span = 0.1,
		color = alpha(dataTypeCols$Ver, 0.8),
		size = 0.2,
		se = FALSE) +
	# Smooth lines for all Types of genes: cross species
	geom_smooth(aes(color = Type, fill = Type), method = "loess", size = 1.5, se = TRUE, span = 0.1) +
	scale_color_manual(values = c(dataTypeCols$HGT, dataTypeCols$Ver), guide = guide_legend(title = "Gene Type")) +
	scale_fill_manual(values = c(dataTypeCols$HGT, dataTypeCols$Ver), guide = FALSE) +
	guides(color = guide_legend(override.aes = list(fill = NA))) +
	lightTheme +
	theme(
		legend.justification = c(1, 1),
		legend.position = c(0.7, 0.825),
		legend.background = element_rect(fill = "white", color = axisCol)
	)


# ------------------------------------------------ #
# Compare Old lHGTs to Recent lHGTS. NB: Recent lHGTS are sparse - hence large SE and inclusion of histogram 
SB_OldvsNew_plot	<- ggplot(SB_OldvsNew_ext, aes(x = BinIndex, y = StrandBias, color = Type)) +
	scale_x_continuous(
		expand = expand_scale(mult = c(0.025, 0.025)),
		name = "Normalized genome position",
		limits = c(0, 1.25),
		breaks = seq(0, 1.25, by = 0.5),
		labels = c("Origin", "Terminus", "Origin"),
		minor_breaks = seq(0, 1.25, by = ((1 / 360) * 30)),
		sec.axis = dup_axis(
			name = NULL,
			breaks = rowMeans(zoneBoundaryList$expandedRange[c("zoneMin", "zoneMax")]),
			labels = zoneBoundaryList$expandedRange$zoneName)) +
	scale_y_continuous(name = "Enrichment for genes on same strand as DnaA", limits = c(-1.5, 1.5)) +
	geom_vline(xintercept = zoneBoundaryList$expandedRange$boundary, color = zoneBoundaryList$boundCol, linetype = "dashed") +	
	geom_rect(data = zoneBoundaryList$expandedRange, aes(xmin = zoneMin, xmax = zoneMax, ymin = -Inf, ymax = Inf), fill = zoneBoundaryList$expandedRange$zoneCol_alpha, inherit.aes = FALSE) +
	# Histogram for Old genes: show per-Bin population (relative)
	geom_histogram(
		data = subset(SB_OldvsNew_ext, Type == "Old" & !(is.na(StrandBias))),
		aes(x = BinIndex,y = (..ncount.. / 4)),
		bins = binNumber,
		fill = alpha(dataTypeCols$Old, 0.4),
		color = alpha(axisCol, 0.8),
		size = 0.2,
		inherit.aes = FALSE) +
	# Histogram for Recent genes: show per-Bin population (relative)
	geom_histogram(
		data = subset(SB_OldvsNew_ext, Type == "Recent" & !(is.na(StrandBias))),
		aes(x = BinIndex,y = -(..ncount.. / 4)),
		bins = binNumber,
		fill = alpha(dataTypeCols$Recent, 0.4),
		color = alpha(axisCol, 0.8),
		size = 0.2,
		inherit.aes = FALSE) +
	# Plot 0-line later to go over bottom of histograms
	geom_hline(yintercept = 0, col = alpha(axisCol, 0.75), size = 1) +
	# Smooth lines for all Types of genes: cross species
	geom_smooth(aes(color = Type, fill = Type), method = "loess", size = 1.5, span = 0.1) +
	scale_color_manual(values = c(dataTypeCols$Old, dataTypeCols$Recent), guide = guide_legend(title = "Gene Type")) +
	scale_fill_manual(values = c(dataTypeCols$Old, dataTypeCols$Recent), guide = FALSE) +
	guides(color = guide_legend(override.aes = list(fill = NA))) +
	lightTheme +
	theme(
		legend.justification = c(1, 1),
		legend.position = c(0.7, 0.825),
		legend.background = element_rect(fill = "white", color = axisCol)
	)



# ------------------------------------------------------------------------------------- #
# Write plots
quartz(width = 18, height = 10)
print(SB_AllGenes_plot)
quartz.save(file = file.path(SBiasFigPath, "All_genes.pdf"), type = "pdf", dpi = 300)
invisible(dev.off())

quartz(width = 18, height = 10)
print(SB_AllvsVer_plot)
quartz.save(file = file.path(SBiasFigPath, "AllvsVertical.pdf"), type = "pdf", dpi = 300)
invisible(dev.off())

quartz(width = 18, height = 10)
print(SB_VervsHGT_plot)
quartz.save(file = file.path(SBiasFigPath, "VerticalvslHGT.pdf"), type = "pdf", dpi = 300)
invisible(dev.off())

quartz(width = 18, height = 10)
print(SB_OldvsNew_plot)
quartz.save(file = file.path(SBiasFigPath, "OldvsRecent_HGT.pdf"), type = "pdf", dpi = 300)
invisible(dev.off())









# ------------------------------------------------ #
# ------------------------------------------------ #
# ------------------------------------------------ #
# Figure 2E

SB_VerOldNew		<- bind_rows(list(byTypeStrandBias_list$Ver, byTypeStrandBias_list$Old, byTypeStrandBias_list$Recent))

SB_VerOldNew_ext	<- SB_VerOldNew %>%
	subset(BinIndex <= 0.125 | BinIndex >= 0.875) %>%
	mutate(BinIndex = case_when(
		BinIndex < 0.5 ~ BinIndex + 1,
		BinIndex > 0.5 ~ BinIndex - 1)) %>%
	bind_rows(SB_VerOldNew, .) %>%
	arrange(BinIndex) %>% 
	mutate(Type = factor(Type, levels = dataTypes_withAge))


# Compare Old lHGTs to Recent lHGTS. NB: Recent lHGTS are sparse - hence large SE and inclusion of histogram 
SB_VerOldNew_plot	<- ggplot(SB_VerOldNew_ext, aes(x = BinIndex, y = StrandBias, color = Type)) +
	scale_x_continuous(
		expand = expand_scale(mult = c(0.01, 0.01)),
		name = "Normalized genome position",
		breaks = seq(0, 1, by = 0.5),
		labels = c("Origin", "Terminus", "Origin"),
		minor_breaks = seq(0, 1, by = ((1 / 360) * 30)),
		sec.axis = dup_axis(
			name = NULL,
			breaks = rowMeans(zoneBoundaryList$fullRange[c("zoneMin", "zoneMax")]),
			labels = zoneBoundaryList$fullRange$zoneName)) +
	scale_y_continuous(
		name = "Enrichment for genes on same strand as DnaA",
		limits = c(-1.5, 1.5)) +
	geom_rect(
		data = zoneBoundaryList$fullRange,
		aes(xmin = zoneMin, xmax = zoneMax, ymin = -Inf, ymax = Inf),
		fill = zoneBoundaryList$fullRange$zoneCol_alpha,
		inherit.aes = FALSE) +
	# Plot 0-line later to go over bottom of histograms
	geom_hline(yintercept = 0, col = alpha(axisCol, 0.75), size = 1) +
	# Smooth lines for all Types of genes: cross species
	geom_smooth(aes(color = Type, fill = Type), method = "loess", size = 1.5, span = 0.1, alpha = 0.2) +
	scale_color_manual(values = c(dataTypeCols$Ver, dataTypeCols$Old, dataTypeCols$Recent), guide = guide_legend(title = "Gene Type")) +
	scale_fill_manual(values = c(dataTypeCols$Ver, dataTypeCols$Old, dataTypeCols$Recent)) +
	guides(color = guide_legend(override.aes = list(fill = NA))) +
	coord_cartesian(xlim = c(0, 1), expand = FALSE) +
	lightTheme +
	theme(
		axis.ticks = element_blank(),
		panel.grid.minor.x = element_blank(),
		legend.justification = c(1, 1),
		legend.position = c(0.7, 0.825),
		legend.background = element_rect(fill = "white", color = axisCol)
	)

quartz(width = 18, height = 10)
print(SB_VerOldNew_plot)
quartz.save(file = file.path(SBiasFigPath, "SB_VerOldNew_plot.pdf"), type = "pdf", dpi = 300)
invisible(dev.off())





