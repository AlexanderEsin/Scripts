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

# Read in subDivision data (bySubdivCOGCompare function requires the colors)
subDivisionKey_data	<- readRDS(file.path(positionData_path, "subDivisionKeyData.rds"))
subDivison_cols		<- subDivisionKey_data$subDivison_cols

message("\rReading in data... done\n")

# ------------------------------------------------------------------------------------- #
# Plotting

# Output path for figures
COGvSpaceFig_path	<- file.path(figureOutput_path, "Space_COGs")

# Quartz plotting options common to this script
quartz.options(canvas = "#333233", bg = "#333233", type = "png", dpi = 300)

# ------------------------------------------------------------------------------------- #
# Compare the COG distributions between the lHGT and Vertical datasets

lHGT_bySubdiv		<- perTypeCOG_data$lHGT$bySubdivision
Ver_bySubdiv		<- perTypeCOG_data$Ver$bySubdivision

# NB subDivision_list is inheritted from masterVariables.R
lHGTvsVer_COGcompare	<- bySubdivCOGCompare(
	bySuvdiv_dataType_A = lHGT_bySubdiv,
	bySuvdiv_dataType_B = Ver_bySubdiv,
	subDivision_list = subDivision_list,
	clusterBy = "A",
	subDivison_cols = subDivison_cols)

# ------------------------------------------------------------------------------------- #
# Compare the COG distributions between the lHGT Old and Recent datasets
Old_bySubdiv		<- perTypeCOG_data$Old$bySubdivision
Recent_bySubdiv		<- perTypeCOG_data$Recent$bySubdivision

OldvsRcent_COGcompare	<- bySubdivCOGCompare(
	bySuvdiv_dataType_A = Old_bySubdiv,
	bySuvdiv_dataType_B = Recent_bySubdiv,
	subDivision_list = subDivision_list,
	clusterBy = "A",
	subDivison_cols = subDivison_cols)


# ------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------- #
# Draw plots

# Plot the lHGT vs Vertical data - clustered by the lHGT distributions
quartz(width = 18, height = 10, file = file.path(COGvSpaceFig_path, "VerticalvslHGT_clusterBylHGT.png"))
pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 1, heights = c(0.5, 1.5))))

pushViewport(viewport(layout.pos.row = 1))
print(lHGTvsVer_COGcompare$ClusterDendro, newpage = FALSE)
popViewport()

pushViewport(viewport(layout.pos.row = 2))
print(lHGTvsVer_COGcompare$ComparisonBarplot, newpage = FALSE)
popViewport()
invisible(dev.off())


# ------------------------------------------------------------------------------------- #

## Plot the lHGT Old vs Recent data - clustered by the Old distributions
quartz(width = 18, height = 10, file = file.path(COGvSpaceFig_path, "OldvsRecent_clusterByOld.png"))
pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 1, heights = c(0.5, 1.5))))

pushViewport(viewport(layout.pos.row = 1))
print(OldvsRcent_COGcompare$ClusterDendro, newpage = FALSE)
popViewport()

pushViewport(viewport(layout.pos.row = 2))
print(OldvsRcent_COGcompare$ComparisonBarplot, newpage = FALSE)
popViewport()
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

