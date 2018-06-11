#!/usr/bin/env Rscript

# Load master variables and HGT position functions
invisible(sapply(HGTPos.all, source, .GlobalEnv))

# Any libraries explicitly used in the script
require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load("plyr", "dplyr", "reshape2", "grid", "ggplot2", "ggrepel", "ggdendro", "ggpubr", "wesanderson", "GGally")


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
# ------------------------------------------------------------------------------------- #
# Important variables
minXferNum		<- 100
contPalette		<- colorRampPalette(wes_palette("Zissou1"))

# ------------------------------------------------------------------------------------- #
# Plotting

# Output path for figures
byBranchFig_path	<- file.path(figureOutput_path, "byBranch")
byZoneFig_path		<- file.path(figureOutput_path, "byZone")
if (!dir.exists(byBranchFig_path)) dir.create(byBranchFig_path)
if (!dir.exists(byZoneFig_path)) dir.create(byZoneFig_path)

# Quartz plotting options common to this script
quartz.options(canvas = "white", bg = "white")

# ------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------- #


# AG trees (from subgroup data)
AG_conTime_tree		<- subgroupData$AG_conTime_tree
AG_binomTime_tree	<- subgroupData$AG_binomTime_tree

# Get a list of unique branches
uniqueBranches			<- unique(perTypeData$lHGT$'3'$allPosData$recepEdge)

# ------------------------------------------------------------------------------------- #

# For each recepEdge (as defined in the lHGT input data), get the transfer data - and circularise the start positions
perBranchData_list	<- lapply(uniqueBranches, function(branch) {
	# All transfered genes at this branch
	perBranchAll	<- perTypeData$lHGT[[hgtPenalty]]$allPosData[which(perTypeData$lHGT[[hgtPenalty]]$allPosData$recepEdge == branch),]
	# Circularise the start positions
	hgtCircStart	<- circular(perBranchAll$relGeneStart * (2 * pi), zero = pi / 2, rotation = "clock", modulo = "2pi")
	return(list(allData = perBranchAll, CircStart = hgtCircStart))
})
names(perBranchData_list)	<- uniqueBranches

## Ignore any branches with less than N transfers
perBranchMinCutoff	<- lapply(perBranchData_list, function(branch) ifelse(length(branch$CircStart) >= minXferNum, return(branch), NA))
perBranchTrim		<- perBranchMinCutoff[!is.na(perBranchMinCutoff)]

## Extract the circular start data
perBranchCircTrim	<- lapply(perBranchTrim, function(branch) return(branch$CircStart))

## Apply the summary.circular_AE function to get, amongst other stats, the Mean and Rho values for each branch
perBranchCircSummary_list	<- lapply(1:length(perBranchCircTrim), function(index) {
	# Get the branch, and summarise
	branchCircSum	<- summary.circular_AE(perBranchCircTrim[[index]])
	# Get the name of the branch
	branchName		<- names(perBranchCircTrim)[index]
	# Return data frame
	BranchCircSum_df	<- cbind(branch = branchName, as.data.frame(t(branchCircSum)), stringsAsFactors = FALSE)
	return(BranchCircSum_df)
})

## Bind the dataframes to make df with all stats per branch
perBranchCircSummary_df	<- bind_rows(perBranchCircSummary_list)

## Add inNode (to match with the tree labels), Index (to make the tree labels), and log(Number of Transfers) to the dataframe
perBranchCircSummary_df$inNode	<- str_split(perBranchCircSummary_df$branch, " ", simplify = TRUE)[,2]
perBranchCircSummary_df$Index	<- rownames(perBranchCircSummary_df)
perBranchCircSummary_df$logN	<- log(perBranchCircSummary_df$n)

## Adjust Index for the root branch to say "Root"
perBranchCircSummary_df$Index[which(perBranchCircSummary_df$inNode == 1375)]	<- "Root"

## Plot the per-Branch compartment population of lHGT genes. These are clustered by distribution
perBranchAvHGT_plot	<- ggplot(perBranchCircSummary_df, aes(x = Mean, y = Rho, color = logN, label = Index)) +
	scale_x_continuous(labels = c("Origin", "Terminus"), breaks = c(0, pi), limits = c(0, 2 * pi)) +
	coord_polar("x") +
	geom_point(size = 1.5) +
	scale_color_gradientn(colours = contPalette(20)) +
	geom_text_repel(aes(label = Index), point.padding = 0.2, min.segment.length = 0.2) +
	lightTheme +
	theme(
		axis.title.x = element_blank(),
	)


# ------------------------------------------------------------------------------------- #
# Prepare dataframes to assign transfers to specific branches on the AG_timeTree

AG_conTime_as4		<- phylo4(AG_conTime_tree)
edgeLabels_df		<- data.frame(E1 = AG_conTime_tree$edge[,1], E2 = AG_conTime_tree$edge[,2], branchLength = AG_conTime_tree$edge.length, inNode = NA, stringsAsFactors = FALSE)
nodeLabels_df		<- data.frame(NodeLabels = labels(AG_conTime_as4, "all"), stringsAsFactors = FALSE)

# Each edge (branch) of the tree gets assigned an "inNode" - the node at which the transfer was predicted to arrive
for (index in 1:nrow(nodeLabels_df)) {
	node	<- nodeLabels_df[index, 1]
	rowID	<- which(nodeLabels_df == node)
	edgeLabels_df$inNode[which(edgeLabels_df$E2 == rowID)]	<- node
}

##We will use the edgeLabels_df to label the tree to match the plot, add the Index labels to match the dataframe
edgeLabels_df$Index		<- NA
for (i in 1:nrow(perBranchCircSummary_df)) {
	row		<- perBranchCircSummary_df[i,]
	edgeLabels_df$Index[which(edgeLabels_df$inNode == row$inNode)]	<- row$Index
}

# Add the Included column to color tree branches by whether they are in the Mean/Rho plot
edgeLabels_df$Included	<- unlist(lapply(edgeLabels_df$Index, function(x) ifelse(is.na(x), 0, 1)))

# ------------------------------------------------------------------------------------- #

perBranchZone_list <- lapply(1:length(perBranchTrim), function(branchData_index) {

	branchData	<- perBranchTrim[[branchData_index]]
	branchName	<- names(perBranchTrim)[branchData_index]
	branchIndex	<- perBranchCircSummary_df$Index[which(perBranchCircSummary_df$branch == branchName)]

	byBranch_zoned		<- splitByZone(data = branchData$allData, zones = zoneBoundaryList$halfGenomeRange, variable = NULL)
	byBranch_zone_df	<- as.data.frame(table(byBranch_zoned$zoneName))

	names(byBranch_zone_df)	<- c("zoneName", "geneNumber")
	byBranch_zone_df$branchName		<- branchName
	byBranch_zone_df$branchIndex	<- branchIndex
	return(list(fullData = byBranch_zoned, table = byBranch_zone_df))
})
perBranchZone_full	<- bind_rows(lapply(perBranchZone_list, function(element) return(element$fullData)))
perBranchZone_df 	<- bind_rows(lapply(perBranchZone_list, function(element) return(element$table)))


## Prepare data to cluster 
clusterRecast_data				<- dcast(perBranchZone_df, branchIndex ~ zoneName, value.var = "geneNumber")
rownames(clusterRecast_data)	<- clusterRecast_data$branchIndex
clusterRecastProp_data			<- sweep(clusterRecast_data[,-1], 1, rowSums(clusterRecast_data[,-1]), "/")
clusterRecastProp_dist			<- dist(as.matrix(clusterRecastProp_data))

clusterZones_dendro	<- dendro_data(hclust(clusterRecastProp_dist, method = "ward.D2"))


## Plot the clustering dendrogram
perBranchZone_cluster	<- ggplot() +
	geom_segment(data = clusterZones_dendro$segments, aes(x = x, y = y, xend = xend, yend = yend), col = axisCol, size = 1) +
	lightTheme +
	theme(
		plot.margin = unit(c(0, 1.45, 0, 1.45), "cm"),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		axis.line = element_blank(),
		axis.text = element_blank(),
		axis.title = element_blank(),
		axis.ticks = element_blank()
	)

# Factor zoneName and branchIndex for plotting.
# For the Index - this means we plot in the same order as the cluster dendrogram
perBranchZone_df$zoneName		<- factor(perBranchZone_df$zoneName, levels = zoneBoundaryList$halfGenomeRange$zoneName)
perBranchZone_df$branchIndex	<- factor(perBranchZone_df$branchIndex, levels = clusterZones_dendro$labels$label)

## Plot the proportional barplots. Each barplot = a single branch, showing proportion of lHGTs in each genome compartment/
perBranchZone_barplot	<- ggplot(data = perBranchZone_df, aes(x = branchIndex, y = geneNumber, fill = zoneName, label = as.character(geneNumber))) +
	geom_bar(stat = "identity", position = "fill") +
	geom_text(position = position_fill(vjust = 0.5), color = "black") +
	scale_fill_manual(values = zoneBoundaryList$halfGenomeRange$zoneColbyName, guide = FALSE) +
	lightTheme +
	theme(
		plot.margin = unit(c(0, 2, 0.5, 2), "cm"),
		panel.grid.major.y = element_blank(),
		panel.grid.major.x = element_line(size = 1, color = axisCol),
		panel.grid.minor.y = element_blank(),
		axis.line = element_blank(),
		axis.text.y = element_blank(),
		axis.title.y = element_blank(),
		axis.ticks = element_blank()
	)

# ------------------------------------------------------------------------------------- #

# # Make a zone key plot

# Zone data
zoneBoundaryList$fullRange$zoneName	<- factor(zoneBoundaryList$fullRange$zoneName, levels = zoneBoundaryList$fullRange$zoneName[1:5])

# Prepare a bare circular plot with the subcompartment zones
zoneKey_plot		<- ggplot(data = zoneBoundaryList$fullRange, aes(xmin = (2 * pi * zoneMin), xmax = (2 * pi * zoneMax), ymin = -Inf, ymax = Inf, fill = zoneName, col = NA)) +
	coord_polar("x") +
	scale_x_continuous(labels = c("Origin", "Terminus"), breaks = c(0, pi), limits = c(0, 2 * pi)) +
	scale_y_continuous(limits = c(0, 1)) +
	geom_rect(col = NA, size = 0) +
	scale_fill_manual(values = zoneBoundaryList$fullRange$zoneColbyName) +
	lightTheme +
	theme(
		panel.grid = element_blank(),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		axis.ticks = element_blank(),
		axis.text.y = element_blank(),
		axis.text.x = element_text(size = 18),
		panel.border = element_blank(),
		legend.key.size = unit(1, "cm"),
		legend.title = element_blank()
	)

# ------------------------------------------------------------------------------------- # 
# ------------------------------------------------------------------------------------- # 
# Write out data plots
# ------------------------------------------------------------------------------------- # 

# Plot the Compartment (Zone) key wheel
quartz(width = 5, height = 4)
print(zoneKey_plot)
quartz.save(file = file.path(byBranchFig_path, "zoneKey.pdf"), type = "pdf", dpi = 100)
quartz.save(file = file.path(byZoneFig_path, "zoneKey.pdf"), type = "pdf", dpi = 100)
invisible(dev.off())

# ------------------------------------------------------------------------------------- # 

# Plot the dendro tree together with the barplots showing lHGT proportion in different genomic compartments
quartz(width = 16, height = 10)
pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 1, heights = c(0.5, 1.5))))

pushViewport(viewport(layout.pos.row = 1))
print(perBranchZone_cluster, newpage = FALSE)
popViewport()

pushViewport(viewport(layout.pos.row = 2))
print(perBranchZone_barplot, newpage = FALSE)
popViewport()

quartz.save(file = file.path(byBranchFig_path, "byBranch_zone.pdf"), type = "pdf", dpi = 100)
invisible(dev.off())

# ------------------------------------------------------------------------------------- # 

# Plot the perBranch HGT mean weight side-by-side with the timeTree labelled with the branches
quartz(width = 9, height = 8)
print(perBranchAvHGT_plot)
quartz.save(file = file.path(byBranchFig_path, "byBranch_circAV.pdf"), type = "pdf", dpi = 100)
invisible(dev.off())

# ------------------------------------------------------------------------------------- # 
# Plot reference tree - colour branches that are used in analysis in red, otherwise grey. Label included branches with the Index number

# Color for the binomial tree
Royal1Pal	<- colorRampPalette(wes_palette("Royal1")[1:2])

quartz(width = 12, height = 8)
plotBranchbyTrait_AE(AG_binomTime_tree, edgeLabels_df$Included, method = "edges", legend = FALSE, palette = Royal1Pal, title = "Included", tip.color = "black", cex = 0.8)
# Annotate tree with Index
for (rowIndex in 1:nrow(edgeLabels_df)) {
	row		<- edgeLabels_df[rowIndex,]
	E2		<- row$E2

	edgeNum	<- which(AG_conTime_tree$edge[,2] == E2)

	if (!is.na(row$Index)) {
		edgelabels(row$Index, edgeNum, adj = c(0.5, -0.25), bg = "#white", frame = "none", col = "black", cex = 1)
	}
}
quartz.save(file = file.path(byBranchFig_path, "indexTree.pdf"), type = "pdf", dpi = 100)
invisible(dev.off())








# ------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------- #


# Boxplot of fraction of HGT into a given category
perBranchZone_forStats	<- dcast(perBranchZone_df, branchIndex ~ zoneName, value.var = "geneNumber")
perBranchZone_forStats	<- cbind(Index = as.character(perBranchZone_forStats$branchIndex), sweep(perBranchZone_forStats[,-1], 1, rowSums(perBranchZone_forStats[,-1]), "/"), stringsAsFactors = FALSE)
perBranchZone_forStats	<- recast(perBranchZone_forStats, variable ~ Index, id.var = "Index")
perBranchZone_forStats	<- melt(perBranchZone_forStats, id.vars = "variable")
names(perBranchZone_forStats)	<- c("Zone", "Index", "Proportion")

# Add branch length to this dataframe
perBranchZone_forStats	<- bind_rows(lapply(unique(perBranchZone_forStats$Index), function(branchIndex) {
	subset_df		<- perBranchZone_forStats[which(perBranchZone_forStats$Index == branchIndex),]
	branchLength	<- ifelse(identical(as.character(branchIndex), "Root"), NA, edgeLabels_df$branchLength[which(edgeLabels_df$Index == branchIndex)]) 
	subset_df$branchLength	<- branchLength
	return(subset_df)
}))

# Factorise the Zones
levels(perBranchZone_forStats$Zone)	<- zoneBoundaryList$halfGenomeRange$zoneName

# Prepare a boxplot showing the relative distribution of lHGT genes across Genome compartments
perBranchByZone_boxplot	<- ggplot(data = perBranchZone_forStats, aes(x = Zone, y = Proportion, fill = Zone)) +
	geom_boxplot(color = axisCol, size = 0.7) +
	scale_y_continuous(name = "Proportion of lHGT per branch", limits = c(0, round_any(max(perBranchZone_forStats$Proportion), 0.1, ceiling))) +
	scale_x_discrete(name = "Genome Zone") +
	scale_fill_manual(values = zoneBoundaryList$fullRange$zoneColbyName, guide = FALSE) +
	stat_compare_means(size = 6, color = axisCol, geom = "label", fill = "white", label.size = 0.3) + 
	lightTheme + 
	theme(
		panel.grid.minor.y = element_blank()
	)


# Correlating lHGT population in each compartments (per branch) against the length of that branch
# Since the root has no length, we remove the root data
perBranchZone_noRoot	<- perBranchZone_forStats[which(perBranchZone_forStats$Index != "Root"),]
perBranchByZone_cor		<- ggplot(data = perBranchZone_noRoot, aes(x = Proportion, y = branchLength, color = Zone, group = Zone)) +
	scale_y_continuous(name = "Branch Length", limits = c(0, round_any(max(perBranchZone_noRoot$branchLength), 0.1, ceiling))) +
	scale_x_continuous(name = "Proportion of lHGT per branch", limits = c(0, round_any(max(perBranchZone_noRoot$Proportion), 0.1, ceiling))) +
	geom_point(size = 2) +
	scale_color_manual(values = zoneBoundaryList$fullRange$zoneColbyName, guide = FALSE) +
	facet_wrap(~Zone) + 
	stat_cor(method = "spearman", geom = "label", fill = "white", label.x.npc = "right", label.y.npc = "top", size = 4) +
	lightTheme +
	theme(
		panel.spacing = unit(2, "lines"),
		strip.background = element_rect(fill = "transparent", color = axisCol, size = 1),
		strip.text = element_text(color = textCol, size = 12)
	)

# ------------------------------------------------------------------------------------- #

## Pairwise comparison of lHGT increase / decrease per compartment against other compartments

# Recast data to have the compartments as individual columns
perBranchZone_forPairwise	<- dcast(perBranchZone_forStats, Index ~ Zone, value.var = "Proportion")
# Max value
maxVal	<- round_any(max(perBranchZone_forStats$Proportion), 0.1, ceiling)

# Plot pairwise correlations
perBranchZone_plot <- ggpairs(data = perBranchZone_forPairwise, columns = 2:ncol(perBranchZone_forPairwise),
	
	upper = list(continuous = function(data, mapping, ...) {
		cor_fun(data = data, mapping = mapping, method = "pearson", ndp = 2, sz = 6, color = axisCol)
	}),

	lower = list(continuous = function(data, mapping, ...) {
		smooth_lm_fun(data = data, mapping = mapping, smooth.colour = alpha(wes_palette("Moonrise1")[3], 0.5), color = axisCol) +
		scale_y_continuous(limits = c(0, maxVal)) +
		scale_x_continuous(limits = c(0, maxVal)) +
		theme(
			panel.grid.major = element_line(size = 0.2, color = alpha(axisCol, alpha = 0.5)),
			panel.grid.minor = element_blank()
		)
	}),
	
	diag = list(continuous = function(data, mapping, ...) {
		ggally_densityDiag(data = data, mapping = mapping, color = axisCol, size = 1.5) +
		scale_x_continuous(limits = c(0, maxVal)) +
		theme(
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank()
		)
	})
)

# Recolour the diagonals to the zone colours
for (i in 1:5) {
	perBranchZone_plot[i,i]$layers[[1]]$aes_params$colour <- zoneBoundaryList$fullRange$zoneColbyName[i]
}

# Apply a variant of darkTheme
zonePairwise_plot	<- perBranchZone_plot + 
	theme(
		plot.background = element_rect(fill = "white", color = NA),
		panel.background = element_rect(fill = "transparent", color = NA),
		axis.ticks = element_line(size = 0.2, color = axisCol),
		plot.title = element_text(size = 16, hjust = 0.5, color = textCol),
		legend.key = element_rect(fill = "transparent", color = NA),
		legend.background = element_rect(fill = "transparent", color = NA),
		legend.text = element_text(size = 14, color = textCol),
		legend.title = element_text(size = 14, color = textCol),
		legend.title.align = 0.5,
		axis.title = element_text(size = 14, color = textCol),
		axis.text = element_text(size = 14, colour = textCol),
		panel.spacing = unit(2, "lines"),
		strip.background = element_rect(fill = "transparent", color = axisCol, size = 1),
		strip.text = element_text(color = textCol, size = 12)
	)


# ------------------------------------------------------------------------------------- #

#  /// Write out the stats plots /// ###

quartz(width = 12, height = 8)
print(perBranchByZone_boxplot)
quartz.save(file = file.path(byZoneFig_path, "byZone_boxplot.pdf"), type = "pdf", dpi = 100)
invisible(dev.off())

quartz(width = 16, height = 10)
print(perBranchByZone_cor)
quartz.save(file = file.path(byZoneFig_path, "byZone_branchVslHGTCor.pdf"), type = "pdf", dpi = 100)
invisible(dev.off())

quartz(width = 16, height = 16)
print(zonePairwise_plot)
quartz.save(file = file.path(byZoneFig_path, "byZone_pairwiseDependence.pdf"), type = "pdf", dpi = 100)
invisible(dev.off())

