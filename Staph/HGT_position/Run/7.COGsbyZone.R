#!/usr/bin/env Rscript
require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load("circular", "dplyr", "grid", "ggplot2", "ggdendro")

# Load master variables and HGT position functions
invisible(sapply(HGTPos.bac, source, .GlobalEnv))

# ------------------------------------------------------------------------------------- #
# Read in data

# message("\nReading in data...", appendLF = FALSE)

# # Per-type COG data
# perTypeCOG_data		<- readRDS(file.path(positionData_path, "AG_perTypeCOGData.rds"))

# # Zone list
# zoneBoundaryList	<- readRDS(file.path(positionData_path, "AG_zoneBoundaries.rds"))

# message("\rReading in data... done\n")

# ------------------------------------------------------------------------------------- #
# Plotting

# Output path for figures
COGvSpaceFig_path	<- file.path(figureOutput_path, "COGs_byZone")
if(!dir.exists(COGvSpaceFig_path)) dir.create(COGvSpaceFig_path)

# Quartz plotting options common to this script
quartz.options(canvas = "white", bg = "white")

# ------------------------------------------------------------------------------------- #
# Compare the COG distributions between the lHGT and Vertical datasets
# lHGT_byZone		<- perTypeCOG_data$lHGT$byZone
# Ver_byZone		<- perTypeCOG_data$Ver$byZone

# NB subDivision_list is inheritted from masterVariables.R
# byZone_dataType_A = lHGT_byZone;
# 	byZone_dataType_B = Ver_byZone;
# 	zones = zoneBoundaryList$fullRange;
# 	clusterBy = "A"

# lHGTvsVer_COGcompare	<- byZoneCOGCompare(
# 	byZone_dataType_A = lHGT_byZone,
# 	byZone_dataType_B = Ver_byZone,
# 	zones = zoneBoundaryList$fullRange,
# 	clusterBy = "A")

# # Re-plot but allow totally free clustering of all COGs
# lHGTvsVer_COGcompare_free	<- byZoneCOGCompare(
# 	byZone_dataType_A = lHGT_byZone,
# 	byZone_dataType_B = Ver_byZone,
# 	zones = zoneBoundaryList$fullRange)

# ------------------------------------------------------------------------------------- #
# Compare the COG distributions between the lHGT Old and Recent datasets
# Old_bySubdiv		<- perTypeCOG_data$Old$byZone
# Recent_bySubdiv		<- perTypeCOG_data$Recent$byZone

# OldvsRcent_COGcompare	<- byZoneCOGCompare(
# 	byZone_dataType_A = Old_bySubdiv,
# 	byZone_dataType_B = Recent_bySubdiv,
# 	zones = zoneBoundaryList$halfGenomeRange,
# 	clusterBy = "A")

# # Re-plot but allow totally free clustering of all COGs
# OldvsRcent_COGcompare_free	<- byZoneCOGCompare(
# 	byZone_dataType_A = Old_bySubdiv,
# 	byZone_dataType_B = Recent_bySubdiv,
# 	zones = zoneBoundaryList$halfGenomeRange)



# ------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------- #
# Draw plots

# Plot the lHGT vs Vertical data - clustered by the lHGT distributions
# quartz(width = 18, height = 10)
# pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 1, heights = c(0.5, 1.5))))

# pushViewport(viewport(layout.pos.row = 1))
# print(lHGTvsVer_COGcompare$ClusterDendro, newpage = FALSE)
# popViewport()

# pushViewport(viewport(layout.pos.row = 2))
# print(lHGTvsVer_COGcompare$ComparisonBarplot, newpage = FALSE)
# popViewport()
# quartz.save(file = file.path(COGvSpaceFig_path, "VerticalvslHGT_clusterBylHGT.pdf"), type = "pdf", dpi = 100)
# invisible(dev.off())

# # Plot the lHGT vs Vertical data - clustered freely
# quartz(width = 20, height = 10)
# pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 1, heights = c(0.5, 1.5))))

# pushViewport(viewport(layout.pos.row = 1))
# print(lHGTvsVer_COGcompare_free$ClusterDendro + theme(plot.margin = unit(c(0, 0.5, 0, 0.5), "cm")), newpage = FALSE)
# popViewport()

# pushViewport(viewport(layout.pos.row = 2))
# print(lHGTvsVer_COGcompare_free$ComparisonBarplot, newpage = FALSE)
# popViewport()
# quartz.save(file = file.path(COGvSpaceFig_path, "VerticalvslHGT_clusterFree.pdf"), type = "pdf", dpi = 100)
# invisible(dev.off())


# ------------------------------------------------------------------------------------- #

# Plot the lHGT Old vs Recent data - clustered by the Old distributions
# quartz(width = 18, height = 10)
# pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 1, heights = c(0.5, 1.5))))

# pushViewport(viewport(layout.pos.row = 1))
# print(OldvsRcent_COGcompare$ClusterDendro, newpage = FALSE)
# popViewport()

# pushViewport(viewport(layout.pos.row = 2))
# print(OldvsRcent_COGcompare$ComparisonBarplot, newpage = FALSE)
# popViewport()
# quartz.save(file = file.path(COGvSpaceFig_path, "OldvsRecent_clusterByOld.pdf"), type = "pdf", dpi = 100)
# invisible(dev.off())


# # Plot the lHGT Old vs Recent data - clustered freely
# quartz(width = 20, height = 10)
# pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 1, heights = c(0.5, 1.5))))

# pushViewport(viewport(layout.pos.row = 1))
# print(OldvsRcent_COGcompare_free$ClusterDendro, newpage = FALSE)
# popViewport()

# pushViewport(viewport(layout.pos.row = 2))
# print(OldvsRcent_COGcompare_free$ComparisonBarplot, newpage = FALSE)
# popViewport()
# quartz.save(file = file.path(COGvSpaceFig_path, "OldvsRecent_clusterFree.pdf"), type = "pdf", dpi = 100)
# invisible(dev.off())

# Interestingly, in the old vs recent, we find that the distribution of 
# metabolic genes (towards origin) - G, C, E (and "-") are not
# significantly different between the age classes.

# ------------------------------------------------------------------------------------- #
# Optional

## Extracting and plotting just the circular weight position plots
# invisible(lapply(seq(1:6), function(plotIndex) {
# 	pushViewport(viewport(layout.pos.col = plotIndex))
# 	print(perTypeCOG_data[[plotIndex]]$plot, newpage = FALSE)
# 	popViewport()
# }))



## ------------------------------------------------------------------------------------- ##
## Heatmap COG distribution ##

COG_minCutOff = 50

lHGT_data	<- perTypeData$lHGT$`4`$allPosData[-which(lapply(perTypeData$lHGT$`4`$allPosData$COGcat, length) > 1),]
lHGT_data$COGcat	<- unlist(lHGT_data$COGcat)

COGsAboveCutOff	<- lHGT_data %>%
	group_by(COGcat) %>%
	summarise(totalGenes = n()) %>%
	subset(totalGenes >= COG_minCutOff)

# Remove those COGs that fall below COG number threshold
lHGT_byZone_cutoff	<- lHGT_data %>%
	subset(COGcat %in% COGsAboveCutOff$COGcat)

# ------------------------------------------------------------------------------------- #
# Calculate the total size of the zones
# zoneBoundary_adj	<- zoneBoundaryList$fullRange %>%
# 	mutate(zoneSize = (zoneMax - zoneMin) * 2) %>%
# 	select(c(zoneName, zoneSize)) %>%
# 	mutate(zoneName = factor(zoneName, levels = levels(lHGT_byZone_cutoff$zone)))

# # Adj the "proportion" of HGTs in the zones by the zone size
# lHGT_byZone_adj		<- lHGT_byZone_cutoff %>%
# 	left_join(zoneBoundary_adj, by = c("zone" = "zoneName")) %>%
# 	mutate(geneDensity = numObsv / zoneSize) %>%
# 	group_by(COGcat) %>%
# 	mutate(zoneExp = zoneSize * sum(geneDensity)) %>%
# 	mutate(scaledProp = (geneDensity / sum(geneDensity)) - (zoneExp / sum(zoneExp)))


# # Cluster data
# densityRecast	<-  dcast(lHGT_byZone_adj, COGcat ~ zone, value.var = "Proportion")
# rownames(densityRecast)	<- densityRecast[,1]
# densityRecast	<- densityRecast[,-1]
# distMat <- dist(as.matrix(densityRecast))
# clusterCompartments_dendro	<- dendro_data(hclust(distMat, method = "ward.D2"))
# ggdendrogram(clusterCompartments_dendro, rotate = T)

# ------------------------------------------------------------------------------------- #
# Prepare COG positions
# filterCOGs	<- lapply(COGsAboveCutOff$COGcat, function(COG) {
# 	COG_subset	<- perTypeCOG_data$lHGT$perCOGdata[[COG]]$allData
# 	unlistCOG	<- COG_subset %>% mutate(COGcat = unlist(COGcat))
# 	return(unlistCOG)
# })
# filterCOG_df	<- bind_rows(filterCOGs)

# filterCOG_df$COGcat	<- factor(filterCOG_df$COGcat, levels = clusterCompartments_dendro$labels$label)

sigmoid = function(x) {
   1 / (1 + exp(-x))
}

COGHeatDistribution_plot	<- ggplot(data = lHGT_byZone_cutoff, aes(x = relGeneStart, y = COGcat, group = COGcat)) +
	facet_wrap(~COGcat, nrow = length(lHGT_byZone_cutoff$COGcat), strip.position = "right", scales = "free_y") +
	stat_density(aes(fill = 10^(..scaled..)), position = "identity", geom = "tile", n = 2000, adjust = 1/4) +
	scale_fill_gradientn(colours = c("white", dataTypeCols$HGT)) +
	# scale_fill_gradientn(colours = c(dataTypeCols$HGT, "white")) +
	# geom_vline(xintercept = zoneBoundaryList$fullRange$boundary, color = boundCol, linetype = "dashed") +
	# Zone colouring
	# geom_rect(data = zoneBoundaryList$fullRange, aes(xmin = zoneMin, xmax = zoneMax, ymin = -Inf, ymax = Inf), fill = alpha(rep(zoneBoundaryList$fullRange$zoneCol_alpha, length(levels(filterCOG_df$COGcat))), 0.1), inherit.aes = FALSE) +
	lightTheme +
	theme(
		panel.grid.major.x = element_blank(),
		panel.grid.major.y = element_blank(),
		panel.grid.minor.x = element_blank(),
		axis.ticks = element_blank(),
		strip.background = element_blank(),
		strip.text.y = element_blank()
	)

quartz(width = 14, height = 6)
print(COGHeatDistribution_plot)
outputFileName	<- file.path(COGvSpaceFig_path, "COG_enrichmentHeatmap.pdf")
invisible(dev.copy2pdf(file = outputFileName))
invisible(dev.off())




COG_G	<- filterCOG_df %>% filter(COGcat == "G")
ggplot(data = COG_G, aes(x = relGeneStart)) +
	# stat_density(aes(y = 1/(1+ exp(-1 * ..density..))), position = "identity", geom = "line", n = 2000, adjust = 1/10) +
	stat_density(aes(y = tanh(..density..)), position = "identity", geom = "line", n = 2000, adjust = 1/10) +
	scale_fill_gradientn(colours = c("white", dataTypeCols$HGT)) +
	geom_vline(xintercept = zoneBoundaryList$fullRange$boundary, color = boundCol, linetype = "dashed") +
	# Zone colouring
	# geom_rect(data = zoneBoundaryList$fullRange, aes(xmin = zoneMin, xmax = zoneMax, ymin = -Inf, ymax = Inf), fill = alpha(rep(zoneBoundaryList$fullRange$zoneCol_alpha, length(levels(filterCOG_df$COGcat))), 0.1), inherit.aes = FALSE) +
	lightTheme +
	theme(
		panel.grid.major.x = element_blank(),
		panel.grid.major.y = element_blank(),
		panel.grid.minor.x = element_blank(),
		axis.ticks = element_blank(),
		strip.background = element_blank(),
		strip.text.y = element_blank()
	)

filterCOG_df %>% group_by(COGcat) %>% summarise(n = n())
