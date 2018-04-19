#!/usr/bin/env Rscript

# Load master variables and HGT position functions
invisible(sapply(HGTPos.all, source, .GlobalEnv))

# Any libraries explicitly used in the script
require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load("dplyr", "reshape2", "grid", "ggplot2", "ggrepel", "ggdendro", "ggpubr", "wesanderson")


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
# ------------------------------------------------------------------------------------- #
# Important variables

lHGTpenalty		<- "4"
minXferNum		<- 100
contPalette		<- colorRampPalette(wes_palette("Zissou1"))

# ------------------------------------------------------------------------------------- #
# Plotting

# Output path for figures
byBranchFig_path	<- file.path(figureOutput_path, "byBranch")
bysubdivFig_path	<- file.path(figureOutput_path, "bySubdiv")

# Quartz plotting options common to this script
quartz.options( canvas = "#333233", bg = "#333233", type = "pdf")

# ------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------- #


# AG trees (from subgroup data)
AG_conTime_tree		<- subgroupData$AG_conTime_tree
AG_binomTime_tree	<- subgroupData$AG_binomTime_tree

# Get a list of unique branches
uniqueBranches			<- unique(perTypeData$lHGT$'3'$allPosData$recepEdge)
numBranchesInTree		<- nrow(AG_conTime_tree$edge)
if (!identical(length(uniqueBranches), numBranchesInTree)) {
	stop("The number of receptor edges in dataset is not equal to number of edges in tree")
}

# ------------------------------------------------------------------------------------- #

# For each recepEdge (as defined in the lHGT input data), get the transfer data - and circularise the start positions
perBranchData_list	<- lapply(uniqueBranches, function(branch) {
	# All transfered genes at this branch
	perBranchAll	<- perTypeData$lHGT[[lHGTpenalty]]$allPosData[which(perTypeData$lHGT[[lHGTpenalty]]$allPosData$recepEdge == branch),]
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
	darkTheme +
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

## We will use the edgeLabels_df to label the tree to match the plot, add the Index labels to match the dataframe
edgeLabels_df$Index		<- NA
for (i in 1:nrow(perBranchCircSummary_df)) {
	row		<- perBranchCircSummary_df[i,]
	edgeLabels_df$Index[which(edgeLabels_df$inNode == row$inNode)]	<- row$Index
}

## Add the Included column to color tree branches by whether they are in the Mean/Rho plot
edgeLabels_df$Included	<- unlist(lapply(edgeLabels_df$Index, function(x) ifelse(is.na(x), 0, 1)))

# Color for the binomial tree
Royal1Pal	<- colorRampPalette(rev(wes_palette("Royal1")[1:2]))


# ------------------------------------------------------------------------------------- #

perBranchSubDivision_list <- lapply(1:length(perBranchTrim), function(branchData_index) {

	branchData		<- perBranchTrim[[branchData_index]]
	branch_name		<- names(perBranchTrim)[branchData_index]
	branch_Index	<- perBranchCircSummary_df$Index[which(perBranchCircSummary_df$branch == branch_name)]

	numSubdivisions			<- length(subDivision_list)
	bySubDivision_branch	<- splitBySubdivision(data = branchData$allData, subDivisions = subDivision_list, variable = NA)
	bySubDivision_df		<- cbind(Branch = rep(branch_name, numSubdivisions), Index = rep(branch_Index, numSubdivisions), bySubDivision_branch, stringsAsFactors = FALSE)
	return(bySubDivision_df)
})
perBranchSubDivision_df <- bind_rows(perBranchSubDivision_list)


## Prepare data to cluster 
clusterRecast_data				<- dcast(perBranchSubDivision_df, Index ~ subDiv, value.var = "numObsv")
rownames(clusterRecast_data)	<- clusterRecast_data$Index
clusterRecastProp_data			<- sweep(clusterRecast_data[,-1], 1, rowSums(clusterRecast_data[,-1]), "/")
clusterRecastProp_dist			<- dist(as.matrix(clusterRecastProp_data))

clusterCompartments_dendro		<- dendro_data(hclust(clusterRecastProp_dist, method = "ward.D2"))


## Plot the clustering dendrogram
perBranchCompartment_cluster	<- ggplot() +
	geom_segment(data = clusterCompartments_dendro$segments, aes(x = x, y = y, xend = xend, yend = yend), col = "#D9D9D9", size = 1) + 
	darkTheme +
	theme(
		plot.margin = unit(c(0, 1.2, 0, 1.2), "cm"),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		axis.text = element_blank(),
		axis.title = element_blank(),
		axis.ticks = element_blank()
	)

# Factorise subDivision and Index for plotting.
# For the Index - this means we plot in the same order as the cluster dendrogram
perBranchSubDivision_df$subDiv	<- factor(perBranchSubDivision_df$subDiv, levels = names(subDivision_list))
perBranchSubDivision_df$Index	<- factor(perBranchSubDivision_df$Index, levels = clusterCompartments_dendro$labels$label)

## Plot the proportional barplots. Each barplot = a single branch, showing proportion of lHGTs in each genome compartment/
perBranchCompartment_barplot	<- ggplot(data = perBranchSubDivision_df, aes(x = Index, y = numObsv, fill = subDiv, label = as.character(numObsv))) +
	geom_bar(stat = "identity", position = "fill") +
	geom_text(position = position_fill(vjust = 0.5), color = "#333233") +
	scale_fill_manual(values = subDivison_cols, guide = FALSE) +
	darkTheme +
	theme(
		plot.margin = unit(c(0, 2, 0.5, 2), "cm"),
		panel.grid.major.y = element_blank(),
		panel.grid.major.x = element_line(size = 1, color = "#D9D9D9"),
		panel.grid.minor.y = element_blank(),
		axis.text.y = element_blank(),
		axis.title.y = element_blank(),
		axis.ticks = element_blank()
	)

# ------------------------------------------------------------------------------------- #

# Make subdivision key plot

# Subdivision data
subDivisionKey_df	<- subDivisionKey_data$subDivisionKey_df
subDivison_cols		<- subDivisionKey_data$subDivison_cols

# Prepare a bare circular plot with the subcompartment zones
subDivisionKey_plot		<- ggplot(data = subDivisionKey_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = col, col = NA)) +
	coord_polar("x") +
	scale_x_continuous(labels = c("Origin", "Terminus"), breaks = c(0, pi), limits = c(0, 2 * pi)) +
	scale_y_continuous(limits = c(0, 1)) +
	geom_rect(col = NA, size = 0) +
	scale_fill_manual(values = subDivison_cols) +
	darkTheme +
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

# Plot the Compartment (Subdivision) key wheel
quartz(width = 5, height = 4, file = file.path(byBranchFig_path, "subDiv_key.pdf"))
print(subDivisionKey_plot)
invisible(dev.off())

# ------------------------------------------------------------------------------------- # 

# Plot the dendro tree together with the barplots showing lHGT proportion in different genomic compartments
quartz(width = 20, height = 14, file = file.path(byBranchFig_path, "byBranch_subdiv.pdf"))
pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 1, heights = c(0.5, 1.5))))

pushViewport(viewport(layout.pos.row = 1))
print(perBranchCompartment_cluster, newpage = FALSE)
popViewport()

pushViewport(viewport(layout.pos.row = 2))
print(perBranchCompartment_barplot, newpage = FALSE)
popViewport()
invisible(dev.off())

# ------------------------------------------------------------------------------------- # 

# Plot the perBranch HGT mean weight side-by-side with the timeTree labelled with the branches
quartz(width = 9, height = 8, file = file.path(byBranchFig_path, "byBranch_circAV.pdf"))
print(perBranchAvHGT_plot)
invisible(dev.off())

# ------------------------------------------------------------------------------------- # 

# Plot reference tree - colour branches that are used in analysis in red, otherwise grey. Label included branches with the Index number
quartz(width = 12, height = 8, file = file.path(byBranchFig_path, "indexTree.pdf"))
plotBranchbyTrait_AE(AG_binomTime_tree, edgeLabels_df$Included, method = "edges", legend = FALSE, palette = Royal1Pal, title = "Included", tip.color = "#D9D9D9", cex = 0.8)
# Annotate tree with Index
for (rowIndex in 1:nrow(edgeLabels_df)) {
	row		<- edgeLabels_df[rowIndex,]
	E2		<- row$E2

	edgeNum	<- which(AG_conTime_tree$edge[,2] == E2)

	if (!is.na(row$Index)) {
		edgelabels(row$Index, edgeNum, adj = c(0.5, -0.25), bg = "#333233", frame = "none", col = "#D9D9D9", cex = 1)
	}
}
invisible(dev.off())













# ------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------- #


# Boxplot of fraction of HGT into a given category
perBranchSubdiv_forStats	<- dcast(perBranchSubDivision_df, Index ~ subDiv, value.var = "numObsv")
perBranchSubdiv_forStats	<- cbind(Index = as.character(perBranchSubdiv_forStats$Index), sweep(perBranchSubdiv_forStats[,-1], 1, rowSums(perBranchSubdiv_forStats[,-1]), "/"), stringsAsFactors = FALSE)
perBranchSubdiv_forStats	<- recast(perBranchSubdiv_forStats, variable ~ Index, id.var = "Index")
perBranchSubdiv_forStats	<- melt(perBranchSubdiv_forStats, id.vars = "variable")
names(perBranchSubdiv_forStats)	<- c("SubDivision", "Index", "Proportion")

# Add branch length to this dataframe
perBranchSubdiv_forStats	<- bind_rows(lapply(unique(perBranchSubdiv_forStats$Index), function(branchIndex) {
	subset_df		<- perBranchSubdiv_forStats[which(perBranchSubdiv_forStats$Index == branchIndex),]
	branchLength	<- ifelse(identical(as.character(branchIndex), "Root"), NA, edgeLabels_df$branchLength[which(edgeLabels_df$Index == branchIndex)]) 
	subset_df$branchLength	<- branchLength
	return(subset_df)
}))

# Factorise the Subdivisions
levels(perBranchSubdiv_forStats$SubDivision)	<- c("Ori", "Near Ori", "Far Ori", "Far Ter", "Near Ter", "Ter")

# Prepare a boxplot showing the relative distribution of lHGT genes across Genome compartments
perBranchBySubdiv_boxplot	<- ggplot(data = perBranchSubdiv_forStats, aes(x = SubDivision, y = Proportion, fill = SubDivision)) +
	geom_boxplot(color = "#D9D9D9") +
	scale_y_continuous(name = "Proportion of lHGT per branch", limits = c(0, 0.5)) +
	scale_x_discrete(name = "Genome Subdivision") +
	scale_fill_manual(values = subDivison_cols, guide = FALSE) +
	stat_compare_means(size = 6, color = "#D9D9D9", geom = "label", fill = "#333233", label.size = 0.3) + 
	darkTheme + 
	theme(panel.grid.minor.y = element_blank())


# Correlating lHGT population in each compartments (per branch) against the length of that branch
# Since the root has no length, we remove the root data
perBranchSubdiv_noRoot		<- perBranchSubdiv_forStats[which(perBranchSubdiv_forStats$Index != "Root"),]
perBranchBySubdiv_cor		<- ggplot(data = perBranchSubdiv_noRoot, aes(x = Proportion, y = branchLength, color = SubDivision, group = SubDivision)) +
	scale_y_continuous(name = "Branch Length", limits = c(0, 0.2)) +
	scale_x_continuous(name = "Proportion of lHGT per branch", limits = c(0, 0.5)) +
	geom_point(size = 2) +
	scale_color_manual(values = subDivison_cols, guide = FALSE) +
	facet_wrap(~SubDivision) + 
	stat_cor(method = "spearman", geom = "text", label.x.npc = "right", label.y.npc = "top", size = 4) +
	darkTheme +
	theme(
		panel.spacing = unit(2, "lines"),
		strip.background = element_rect(fill = "transparent", color = "#D9D9D9", size = 1),
		strip.text = element_text(color = "#D9D9D9", size = 12)
		)

# ------------------------------------------------------------------------------------- #

## Pairwise comparison of lHGT increase / decrease per compartment against other compartments

# Recast data to have the compartments as individual columns
perBranchSubdiv_forPairwise	<- dcast(perBranchSubdiv_forStats, Index ~ SubDivision, value.var = "Proportion")

# Plot pairwise correlations
SubDivisionPairwise_plot <- ggpairs(data = perBranchSubdiv_forPairwise, columns = 2:ncol(perBranchSubdiv_forPairwise),
	
	upper = list(continuous = function(data, mapping, ...) {
		cor_fun(data = data, mapping = mapping, method = "pearson", ndp = 2, sz = 6, color = "#D9D9D9")
	}),

	lower = list(continuous = function(data, mapping, ...) {
		smooth_lm_fun(data = data, mapping = mapping, smooth.colour = "green", color = "#D9D9D9") +
		scale_y_continuous(limits = c(0, 0.5)) +
		scale_x_continuous(limits = c(0, 0.5)) +
		theme(
			panel.grid.major = element_line(size = 0.2, color = alpha("#D9D9D9", alpha = 0.5)),
			panel.grid.minor = element_blank()
		)
	}),
	
	diag = list(continuous = function(data, mapping, ...) {
		ggally_densityDiag(data = data, mapping = mapping, color = "#D9D9D9") +
		scale_x_continuous(limits = c(0, 0.5)) +
		theme(
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank()
		)
	})
)

# Recolour the diagonals to the subdivision colours
for (i in 1:6) {
	SubDivisionPairwise_plot[i,i]$layers[[1]]$aes_params$colour <- subDivison_cols[i]
}

# Apply a variant of darkTheme
SubDivisionPairwise_plot	<- SubDivisionPairwise_plot + 
	theme(
		plot.background = element_rect(fill = "#333233", color = NA),
		panel.background = element_rect(fill = "transparent", color = NA),
		axis.ticks = element_line(size = 0.2, color = "#D9D9D9"),
		plot.title = element_text(size = 16, hjust = 0.5, color = "#D9D9D9"),
		axis.title = element_text(size = 14, color = "#D9D9D9"),
		axis.text = element_text(size = 14, colour = "#D9D9D9"),
		panel.spacing = unit(2, "lines"),
		strip.background = element_rect(fill = "transparent", color = "#D9D9D9", size = 1),
		strip.text = element_text(color = "#D9D9D9", size = 14)
	)


# ------------------------------------------------------------------------------------- #

#  /// Write out the stats plots /// ###

quartz(width = 12, height = 8, file = file.path(bysubdivFig_path, "bySubdiv_boxplot.pdf"))
print(perBranchBySubdiv_boxplot)
invisible(dev.off())

quartz(width = 16, height = 10, file = file.path(bysubdivFig_path, "bySubdiv_branchVslHGTCor.pdf"))
print(perBranchBySubdiv_cor)
invisible(dev.off())

quartz(width = 16, height = 16, file = file.path(bysubdivFig_path, "bySubdiv_pairwiseCor.pdf"))
print(SubDivisionPairwise_plot)
invisible(dev.off())

