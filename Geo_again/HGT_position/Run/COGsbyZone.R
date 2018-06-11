#!/usr/bin/env Rscript
require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load("circular", "dplyr", "grid", "ggplot2", "ggdendro")

# Load master variables and HGT position functions
invisible(sapply(HGTPos.all, source, .GlobalEnv))

# ------------------------------------------------------------------------------------- #
# Read in data

message("\nReading in data...", appendLF = FALSE)

# Per-type COG data
perTypeCOG_data		<- readRDS(file.path(positionData_path, "AG_perTypeCOGData.rds"))

# Zone list
zoneBoundaryList	<- readRDS(file.path(positionData_path, "AG_zoneBoundaries.rds"))

message("\rReading in data... done\n")

# ------------------------------------------------------------------------------------- #
# Plotting

# Output path for figures
COGvSpaceFig_path	<- file.path(figureOutput_path, "COGs_byZone")
if(!dir.exists(COGvSpaceFig_path)) dir.create(COGvSpaceFig_path)

# Quartz plotting options common to this script
quartz.options(canvas = "white", bg = "white")

# ------------------------------------------------------------------------------------- #
# Compare the COG distributions between the lHGT and Vertical datasets
lHGT_byZone		<- perTypeCOG_data$lHGT$byZone
Ver_byZone		<- perTypeCOG_data$Ver$byZone

# NB subDivision_list is inheritted from masterVariables.R
lHGTvsVer_COGcompare	<- byZoneCOGCompare(
	byZone_dataType_A = lHGT_byZone,
	byZone_dataType_B = Ver_byZone,
	zones = zoneBoundaryList$halfGenomeRange,
	clusterBy = "A")

# Re-plot but allow totally free clustering of all COGs
lHGTvsVer_COGcompare_free	<- byZoneCOGCompare(
	byZone_dataType_A = lHGT_byZone,
	byZone_dataType_B = Ver_byZone,
	zones = zoneBoundaryList$halfGenomeRange)

# ------------------------------------------------------------------------------------- #
# Compare the COG distributions between the lHGT Old and Recent datasets
Old_bySubdiv		<- perTypeCOG_data$Old$byZone
Recent_bySubdiv		<- perTypeCOG_data$Recent$byZone

OldvsRcent_COGcompare	<- byZoneCOGCompare(
	byZone_dataType_A = Old_bySubdiv,
	byZone_dataType_B = Recent_bySubdiv,
	zones = zoneBoundaryList$halfGenomeRange,
	clusterBy = "A")

# Re-plot but allow totally free clustering of all COGs
OldvsRcent_COGcompare_free	<- byZoneCOGCompare(
	byZone_dataType_A = Old_bySubdiv,
	byZone_dataType_B = Recent_bySubdiv,
	zones = zoneBoundaryList$halfGenomeRange)



# ------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------- #
# Draw plots

# Plot the lHGT vs Vertical data - clustered by the lHGT distributions
quartz(width = 18, height = 10)
pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 1, heights = c(0.5, 1.5))))

pushViewport(viewport(layout.pos.row = 1))
print(lHGTvsVer_COGcompare$ClusterDendro, newpage = FALSE)
popViewport()

pushViewport(viewport(layout.pos.row = 2))
print(lHGTvsVer_COGcompare$ComparisonBarplot, newpage = FALSE)
popViewport()
quartz.save(file = file.path(COGvSpaceFig_path, "VerticalvslHGT_clusterBylHGT.pdf"), type = "pdf", dpi = 100)
invisible(dev.off())

# Plot the lHGT vs Vertical data - clustered freely
quartz(width = 20, height = 10)
pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 1, heights = c(0.5, 1.5))))

pushViewport(viewport(layout.pos.row = 1))
print(lHGTvsVer_COGcompare_free$ClusterDendro + theme(plot.margin = unit(c(0, 0.5, 0, 0.5), "cm")), newpage = FALSE)
popViewport()

pushViewport(viewport(layout.pos.row = 2))
print(lHGTvsVer_COGcompare_free$ComparisonBarplot, newpage = FALSE)
popViewport()
quartz.save(file = file.path(COGvSpaceFig_path, "VerticalvslHGT_clusterFree.pdf"), type = "pdf", dpi = 100)
invisible(dev.off())


# ------------------------------------------------------------------------------------- #

# Plot the lHGT Old vs Recent data - clustered by the Old distributions
quartz(width = 18, height = 10)
pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 1, heights = c(0.5, 1.5))))

pushViewport(viewport(layout.pos.row = 1))
print(OldvsRcent_COGcompare$ClusterDendro, newpage = FALSE)
popViewport()

pushViewport(viewport(layout.pos.row = 2))
print(OldvsRcent_COGcompare$ComparisonBarplot, newpage = FALSE)
popViewport()
quartz.save(file = file.path(COGvSpaceFig_path, "OldvsRecent_clusterByOld.pdf"), type = "pdf", dpi = 100)
invisible(dev.off())


# Plot the lHGT Old vs Recent data - clustered freely
quartz(width = 20, height = 10)
pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 1, heights = c(0.5, 1.5))))

pushViewport(viewport(layout.pos.row = 1))
print(OldvsRcent_COGcompare_free$ClusterDendro, newpage = FALSE)
popViewport()

pushViewport(viewport(layout.pos.row = 2))
print(OldvsRcent_COGcompare_free$ComparisonBarplot, newpage = FALSE)
popViewport()
quartz.save(file = file.path(COGvSpaceFig_path, "OldvsRecent_clusterFree.pdf"), type = "pdf", dpi = 100)
invisible(dev.off())

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

