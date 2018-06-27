#!/usr/bin/env Rscript

# Load master variables and HGT position functions
invisible(sapply(HGTPos.all, source, .GlobalEnv))

require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load("ape", "phylobase", "phangorn", "geiger", "stringr", "wesanderson", "tidyverse", "Hmisc", "broom")

# ------------------------------------------------------------------------------------- #
# Read in data

refinedHGTs_path	<- file.path(mowgli_path, "Mowgli_output", "Cleaned_events", "HGT_events")

# ------------------------------------------------------------------------------------- #?


# Define and read in the list of all AG tips as they appear in the Mowgli species tree
AG_mowCladoTips_file	<- file.path(mowgli_path, "Inside_group", "AnoxyGeo_mowTips.txt")
AG_mowCladoTips_char	<- read.table(file = AG_mowCladoTips_file, sep = "\n", stringsAsFactors = FALSE)$V1

# Define and read in the taxid <-> binomial translation table
taxidBinom_trans_file	<- file.path(genome_path, "Genome_lists", "Taxid_refinedBinomial_table.tsv")
taxidBinom_trans_df		<- read.table(file = taxidBinom_trans_file, sep = "\t", header = TRUE,stringsAsFactors = FALSE)


# ------------------------------------------------ #

# Define and read in the mowgli species tree
mowSpeciesTree_file	<- file.path(mowgli_path, "Mowgli_output", "Output_3", "1", "outputSpeciesTree.mpr")
if (!file.exists(mowSpeciesTree_file)) stop("Cannot find Mowgli species tree file, try a different file")
mowSpecies_tree		<- read.tree(mowSpeciesTree_file)

# Define and read in the consensus time tree (same topology as mowgli tree)
consensusTree_file	<- file.path(master_path, "Consensus_groups", "Bacillaceae", "Final_trees", "Taxid_labelled", "RAxML_bipartitions.super_tree-50.txt")
consensusTree_tree	<- read.tree(file = consensusTree_file)


# ------------------------------------------------ #
# Using the AG mowgli tips labels extract the AG subtree
AG_mowClado_tree		<- drop.tip(mowSpecies_tree, setdiff(mowSpecies_tree$tip.label, AG_mowCladoTips_char))

# Produce a df with corresponding taxid - mowgli tip entries
taxidMowExtend_df			<- data.frame(
	Full = AG_mowClado_tree$tip.label, 
	str_split(AG_mowClado_tree$tip.label, "_", simplify = TRUE),
	stringsAsFactors = FALSE
)

# Final column is the binomials
taxidMowExtend_df[,4]		<- apply(taxidMowExtend_df, 1, function(row) {
	taxid	<- row[2]
	binom	<- taxidBinom_trans_df$Binomial[which(taxidBinom_trans_df$Taxid == taxid)]
})
# Rename columns 
names(taxidMowExtend_df)	<- c("Full", "Taxid", "Extension", "Binomial")

## Write out the table for future use
write.table(taxidMowExtend_df, file = file.path(timeData_path, "Taxid2MowTip_table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Translate the mowgli clado tree tips to binomials
AG_mowCladoBin_tree				<- AG_mowClado_tree
AG_mowCladoBin_tree$tip.label	<- unlist(lapply(AG_mowCladoBin_tree$tip.label, function(tip) {
	newTipLabel_char	<- taxidMowExtend_df$Binomial[which(taxidMowExtend_df$Full == tip)]
}))

# Translate the mowgli clado tree tips to the mowgli node IDs
AG_mowClado_tree$tip.label	<- unlist(lapply(AG_mowClado_tree$tip.label, function(tip) {
	newTipLabel_char	<- taxidMowExtend_df$Extension[which(taxidMowExtend_df$Full == tip)]
}))

# Transform to phylo4 object
AG_mowClado_as4	<- phylo4(AG_mowClado_tree)


# ------------------------------------------------ #
# Extract the AG subtree
AG_conTime_tree		<- drop.tip(consensusTree_tree, setdiff(consensusTree_tree$tip.label, taxidMowExtend_df$Taxid))

# Relabel the time tree to binomials
AG_conTimeBin_tree	<- AG_conTime_tree
AG_conTimeBin_tree$tip.label	<- unlist(lapply(AG_conTimeBin_tree$tip.label, function(tip) {
	newTipLabel_char		<- taxidMowExtend_df$Binomial[which(taxidMowExtend_df$Taxid == tip)]
}))

# Relabel the time tree to the Mowgli extension
AG_conTime_tree$tip.label	<- unlist(lapply(AG_conTime_tree$tip.label, function(tip) {
	newTipLabel_char		<- taxidMowExtend_df$Extension[which(taxidMowExtend_df$Taxid == tip)]
}))

# Remove the BS values and conver the phylo4
AG_conTime_tree$node.label	<- NULL
AG_conTime_as4		<- phylo4(AG_conTime_tree)


# Double check the clado (Astral) and time (consensus) topologies are the same
cladoVsTime_rfDist	<- RF.dist(AG_mowClado_tree, AG_conTime_tree)


# ------------------------------------------------ #
# Relabel internal timeTree nodes to match the Mowgli labels in the clado tree
AG_conTime_intNodes		<- nodeId(AG_conTime_as4, "internal")

# Get the node labels as a data frame
AG_mowClado_intNodes	<- nodeId(AG_mowClado_as4, "internal")

# For each internal node in the time tree find a corresponding internal clado
# node with the same underlying tips, then relabel that time tree node with the
# corresponding clado (mowgli) label
for (timeIntNode in AG_conTime_intNodes) {
	matched	<- FALSE
	for (cladoIntNode in AG_mowClado_intNodes) {
		if (setequal(tips(AG_conTime_tree, timeIntNode), tips(AG_mowClado_tree, cladoIntNode)) == TRUE) {

			## Get the label of the node from the species-cladogram tree ##
			node_label	<- AG_mowClado_as4@label[cladoIntNode]

			## Index of the internal node label in the time-resolved tree ##
			node_index	<- which(names(nodeLabels(AG_conTime_as4)) == timeIntNode)
			
			## Assign the name of the node to the time-resolved tree ##
			nodeLabels(AG_conTime_as4)[node_index] <- node_label

			cat(paste0("Time-resolved tree node: ", timeIntNode, " maps to species tree node: ", cladoIntNode, "\n"))
			matched	<- TRUE
			break
		}
	}
	if (matched == FALSE) {
		stop(paste0("Time tree node: ", timeIntNode, " did not map to any species tree node\n"))
	}
}

## Write out the time tree (for use in other scripts - e.g. HGT_position)
write.tree(as(AG_conTime_as4, "phylo"), file = file.path(timeData_path, "AG_conTimeTree.tree"))


# ------------------------------------------------ #
# A list of two dfs - the lists of all nodes found in the Clado and Time trees
# These will be matched to the receptor nodes of the HGTs below, and is why we 
# needed to relabel the tips to the Mowgli node IDs - to match transfer nodes.
allNodeLabels_dfs	<- list(
	clado = data.frame(NodeLabels = labels(AG_mowClado_as4, "all"), stringsAsFactors = FALSE),
	time = data.frame(NodeLabels = labels(AG_conTime_as4, "all"), stringsAsFactors = FALSE)
)

# Prepare tree-edge based dfs. Each edge, flanked by two nodes (E1 & E2) will
# have transfers assigned. In the HGT prediction the "receptor node" corresponds
# to the bottom of the edge (E2). The list of node (above) and edges appear in
# an identical order within each tree, so we can assign the transfers to the edges
# using the mowgli node IDs.
allEdgeFinal_dfs	<- list(clado = data.frame(AG_mowClado_tree$edge, Subgroup = NA), time = data.frame(AG_conTime_tree$edge, Subgroup = NA))
allEdgeFinal_dfs	<- lapply(allEdgeFinal_dfs, setNames, c("E1", "E2", "Subgroup"))

# Create columns for each penalty tested (minimal penalty for consistent HGT)
# is T4 (i.e. consistent over penalties 3 and 4)
allPenNames_char	<- lapply(penalty_list, function(penalty) paste0("T", penalty))
allEdgeFinal_dfs	<- lapply(allEdgeFinal_dfs, function(df) {
	for (penName in allPenNames_char) {0 -> df[eval(penName)]}
	return(df)
})

## Process lHGT data
# For each penalty, process the lHGTs into each node, and assign them to an edge
# in the tree-edge based dfs
for (penalty in penalty_list) {

	# Define and read in the consistent lHGT file
	lHGT_name_file	<- paste0("T", penalty, "_full_lHGT_events.tsv")
	lHGT_data_df	<- read.table(file = file.path(refinedHGTs_path, lHGT_name_file), header = TRUE, sep = "\t")

	# Define column to which HGTs are assigned
	colName			<- paste0("T", penalty)
	# Counter for transfers at the AG root (not plotted on tree)
	at_root			<- 0

	# Process the lHGT results line-by-line
	for (rowIndex in 1:nrow(lHGT_data_df)) {

		# Define line and extract the "receptor node"
		row 		<- lHGT_data_df[rowIndex,]
		recNodes	<- str_split(row$Receptor_nodes, pattern = " ", simplify = T)

		# Check whether transfer is into the root of AG
		if (recNodes[2] == 1375) {
			at_root	<- at_root + 1
		# This would be a transfer above the root of AG and shouldn't happen
		} else if (recNodes[2] == 1531) {
			message("Should not get this")
		} else {
			# Pull out row number (indication of phylogeny edge) corresponding to the Mowgli tree label.
			edgeForClado	<- grep(recNodes[2], allNodeLabels_dfs$clado$NodeLabels)
			edgeForTime		<- grep(recNodes[2], allNodeLabels_dfs$time$NodeLabels)

			# Shortcut names for the corresponding edge + penalty values (ongoing count)
			cladoEntry		<- allEdgeFinal_dfs$clado[[colName]][which(allEdgeFinal_dfs$clado$E2 == edgeForClado)]
			timeEntry		<- allEdgeFinal_dfs$time[[colName]][which(allEdgeFinal_dfs$time$E2 == edgeForTime)]

			# Update both Time and Clado trees with an extra transfer at the edge + penalty
			allEdgeFinal_dfs$clado[[colName]][which(allEdgeFinal_dfs$clado$E2 == edgeForClado)]	<-  cladoEntry + 1      	
			allEdgeFinal_dfs$time[[colName]][which(allEdgeFinal_dfs$time$E2 == edgeForTime)]	<-  timeEntry + 1
		}
	}
	message(paste0("Penalty: " , penalty, " - at root transfers: ", at_root))
}


# ------------------------------------------------ #
# Read in subdivisions

# Read in the table constructed based on the publication #
subspeciesGroup_file	<- file.path(genome_path, "Genome_lists", "AG_subspeciesGroups.txt")
subspeciesGroup_data	<- read.table(file = subspeciesGroup_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Get those tips that have a fellow tip in the same subgroup. Isolate all the subgroups that have more than one member #
subspecOnly_data		<- subset(subspeciesGroup_data, duplicated(subspeciesGroup_data$Group) | duplicated(subspeciesGroup_data$Group, fromLast = TRUE))
uniqueSubGroup_char		<- unique(subspecOnly_data$Group)

# Label branches that connect tips belonging to the same subgroup:
# For each group with more than one member, identify the tips. Find all the
# descendant nodes from the common ancestor of the subgroup (they are
# monophyletic). Using those nodes, identify all the edges (branches) that we
# will classify as being a subgroup branch (1). All other branches are separate
# groups (0)
for (subGroup in uniqueSubGroup_char) {

	taxidsInSubGroup	<- as.vector(subspecOnly_data$Taxid[which(subspecOnly_data$Group == subGroup)])
	taxidsToMowTips		<- taxidMowExtend_df$Extension[which(taxidMowExtend_df$Taxid %in% taxidsInSubGroup)]
	allSubGroupNodes	<- as.vector(descendants(AG_conTime_as4, MRCA(AG_conTime_as4, taxidsToMowTips), "all"))

	for (node in allSubGroupNodes) {
		allEdgeFinal_dfs$time$Subgroup[which(allEdgeFinal_dfs$time$E2 == node)]		<- 1
		allEdgeFinal_dfs$clado$Subgroup[which(allEdgeFinal_dfs$clado$E2 == node)]	<- 1
	}
}

# All other edges labelled as 0 (not subgroup)
allEdgeFinal_dfs	<- lapply(allEdgeFinal_dfs, function(df) {
	df$Subgroup[is.na(df$Subgroup)]	<- 0
	return(df)
})

# ------------------------------------------------ #

# Plot a phylogram colored according to the OLD / RECENT colors with a scale bar
quartz(width = 20, height = 10)

OldRecent_pal	<- colorRampPalette(c(dataTypeCols$Old, dataTypeCols$Recent))
plotBranchbyTrait_AE(AG_conTimeBin_tree, allEdgeFinal_dfs$time$Subgroup, method = "edges", palette = OldRecent_pal, title = "Subgroup branches")
add.scale.bar()

quartz.save(file = file.path(timeOutput_path, "AG_25Genome_phyloTree_oldRecent.pdf"), type = "pdf", dpi = 300)
invisible(dev.off())



# ------------------------------------------------ #
# Calculate fraction of all HGT
allEdgeFinal_fract_dfs	<- lapply(allEdgeFinal_dfs, function(df) cbind(df[,1:3], sweep(df[,-1:-3] , 2, colSums(df[,-1:-3]), "/") * 100))
allEdgeFinal_fract_dfs	<- lapply(allEdgeFinal_fract_dfs, function(df) cbind(df, Mean = rowMeans(df[,-1:-3])))


# HGT-per-branch - Clado

# Color palette for HGT / branch 
conPal <- colorRampPalette(wes_palette("Zissou1"))
# Plot mean fraction of transfer per branch on the Clado tree
plotBranchbyTrait_AE(AG_mowCladoBin_tree, allEdgeFinal_fract_dfs$clado$Mean, method = "edges", legend = 5, palette = conPal, title = "Fraction of transfers")
# Anotate the branches with the numeric values
for (index in 1:nrow(allEdgeFinal_fract_dfs$clado)) {
	entry <- allEdgeFinal_fract_dfs$clado[index,]
	AnnotateBranch(entry, index, "Mean")
}


## HGT-per-branch - Time
# Plot mean fraction of transfer per branch on the Time tree
plotBranchbyTrait_AE(AG_conTimeBin_tree, allEdgeFinal_fract_dfs$time$Mean, method = "edges", legend = .1, palette = conPal, title = "Fraction of transfers")

# Anotate the branches with the numeric values
for (index in 1:nrow(allEdgeFinal_fract_dfs$time)) {
	entry <- allEdgeFinal_fract_dfs$time[index,]
	AnnotateBranch(entry, index, "Mean")
}


# ------------------------------------------------ #
# Continuous HGT

## HGT vs Branch Length
allEdgeFinal_dfs$time	<- data.frame(
	allEdgeFinal_dfs$time,
	Mean = allEdgeFinal_fract_dfs$time$Mean,
	BranchLen = AG_conTime_tree$edge.length,
	Index = rownames(allEdgeFinal_dfs$time),
	SubGroupLogic = as.logical(allEdgeFinal_dfs$time$Subgroup),
	stringsAsFactors = FALSE
)


branchVsHGT_plot	<- ggplot(data = allEdgeFinal_dfs$time, mapping = aes(x = BranchLen, y = T4, label = Index, color = SubGroupLogic, group = SubGroupLogic)) +
	geom_point(size = 4) +
	scale_x_continuous(
		name = "Branch length",
		limits = c(0, 0.2)) + 
	scale_y_continuous(
		name = "Number of HGT events per branch",
		limits = c(0, 160)) +
	annotate("label", x = 0.15, y = 120, label = paste("Total HGT events = ", sum(allEdgeFinal_dfs$time$T4))) +
	# geom_smooth(method = 'lm', formula = y~x, se = FALSE) +
	scale_color_manual(values = c(dataTypeCols$Old, dataTypeCols$Recent)) +
	lightTheme +
	theme(
		axis.ticks = element_blank()
	)


branchVsHGT_branchLabel_plot	<- branchVsHGT_plot + geom_text(aes(label = Index), hjust = 0, vjust = -1)


quartz(width = 15, height = 10)
print(branchVsHGT_plot)
quartz.save(file = file.path(timeOutput_path, "AG_branchVsHGT_scatter.pdf"), type = "pdf", dpi = 300)
invisible(dev.off())

quartz(width = 15, height = 10)
print(branchVsHGT_branchLabel_plot)
quartz.save(file = file.path(timeOutput_path, "AG_branchVsHGT_branchLabel_scatter.pdf"), type = "pdf", dpi = 300)
invisible(dev.off())


quartz(width = 20, height = 10)
plotBranchbyTrait_AE(AG_mowCladoBin_tree, allEdgeFinal_dfs$time$Subgroup, method = "edges", palette = OldRecent_pal, title = "Subgroup branches")
for (i in as.vector(allEdgeFinal_dfs$time$Index)) {
	edgelabels(i, as.numeric(i), adj = c(0.5, -0.25), bg = "white", frame = "none", cex = 0.8)
}
quartz.save(file = file.path(timeOutput_path, "AG_25Genome_cladoTree_oldRecent_branchLabel.pdf"), type = "pdf", dpi = 300)
invisible(dev.off())


# ------------------------------------------------ #


perPenaltyCor_list	<- lapply(c(penalty_list, "Mean"), function(penalty) {
	
	if (!identical(penalty, "Mean")) {
		penColName	<- paste0("T", penalty)
	} else {
		penColName	<- penalty
	}
	time_df			<- allEdgeFinal_dfs$time

	# Overall pearson + rsquared
	all_cor			<- tidy(rcorr(time_df$BranchLen, time_df[[penColName]]))
	all_corVal		<- all_cor$estimate
	all_pVal		<- all_cor$p.value
	linear_mod		<- summary(lm(time_df$BranchLen ~ time_df[[penColName]]))
	all_rSqr_val	<- linear_mod$r.squared
	all_rSqr_pval	<- tidy(linear_mod)$p.value[2]

	# Pearson + rsquared (without two obvious outliers)
	outlierIndex	<- as.character(c(1, 8))
	timeNoOut_df	<- time_df[-which(time_df$Index %in% outlierIndex),]

	noOut_cor		<- tidy(rcorr(timeNoOut_df$BranchLen, timeNoOut_df[[penColName]]))
	noOut_corVal	<- noOut_cor$estimate
	noOut_pVal		<- noOut_cor$p.value
	linear_mod		<- summary(lm(timeNoOut_df$BranchLen ~ timeNoOut_df[[penColName]]))
	noOut_rSqr_val	<- linear_mod$r.squared
	noOut_rSqr_pval	<- tidy(linear_mod)$p.value[2]

	# Correlation when subgroup removed
	grpTime_df	<- time_df[which(time_df$SubGroupLogic == FALSE),]

	grp_cor		<- tidy(rcorr(grpTime_df$BranchLen, grpTime_df[[penColName]]))
	grp_corVal	<- grp_cor$estimate
	grp_pVal	<- grp_cor$p.value
	linear_mod		<- summary(lm(grpTime_df$BranchLen ~ grpTime_df[[penColName]]))
	grp_rSqr_val	<- linear_mod$r.squared
	grp_rSqr_pval	<- tidy(linear_mod)$p.value[2]

	# Pearson + rsquared (without two obvious outliers)
	grpTimeNoOut_df	<- grpTime_df[-which(grpTime_df$Index %in% outlierIndex),]

	grpNoOut_Pcor	<- tidy(rcorr(grpTimeNoOut_df$BranchLen, grpTimeNoOut_df[[penColName]]))
	grpNoOut_corVal	<- grpNoOut_Pcor$estimate
	grpNoOut_pVal	<- grpNoOut_Pcor$p.value
	linear_mod		<- summary(lm(grpTimeNoOut_df$BranchLen ~ grpTimeNoOut_df[[penColName]]))
	grpNoOut_rSqr_val	<- linear_mod$r.squared
	grpNoOut_rSqr_pval	<- tidy(linear_mod)$p.value[2]


	outputPearson	<- data.frame(
		Penalty = penalty,
		Test = "Pearson",
		all = all_corVal,
		all_pval = all_pVal,
		noOut = noOut_corVal,
		noOut_pval = noOut_pVal,
		group = grp_corVal,
		group_pval = grp_pVal,
		groupNoOut = grpNoOut_corVal,
		groupNoOut_pval = grpNoOut_pVal,
		stringsAsFactors = FALSE)

	outputRsqrd		<- data.frame(
		Penalty = penalty,
		Test = "Rsquared",
		all = all_rSqr_val,
		all_pval = all_rSqr_pval,
		noOut = noOut_rSqr_val,
		noOut_pval = noOut_rSqr_pval,
		group = grp_rSqr_val,
		group_pval = grp_rSqr_pval,
		groupNoOut = grpNoOut_rSqr_val,
		groupNoOut_pval = grpNoOut_rSqr_pval,
		stringsAsFactors = FALSE)

	return(bind_rows(list(outputPearson, outputRsqrd)))
})

perPenaltyCor_df	<- bind_rows(perPenaltyCor_list)
perPenaltyCor_df	<- perPenaltyCor_df %>% mutate_if(is.numeric, funs(signif(., digits = 3)))


perPenaltyCor_melt	<- melt(perPenaltyCor_df, id.vars = c("Penalty", "Test"), variable.name = "Subset", value.name = "StatisticValue")
perPenaltyCor_melt$Penalty	<- factor(perPenaltyCor_melt$Penalty, levels = unique(perPenaltyCor_melt$Penalty))
# Filter out the pvalues for plotting
perPenaltyCor_melt 	<- perPenaltyCor_melt %>% filter(!str_detect(Subset, "pval"))

perPenaltyCor_plot	<- ggplot(data = perPenaltyCor_melt, mapping = aes(x = Penalty, y = StatisticValue, color = Subset, group = Subset)) +
	geom_point(size = 3) +
	geom_line(size = 2) +
	scale_y_continuous(limits = c(0, 1)) +
	facet_wrap(~Test) +
	scale_color_manual(values = rev(wes_palette("BottleRocket2")[1:4])) +
	theme(panel.background = element_blank(), plot.title = element_text(hjust = 0.5), axis.line = element_line(colour = "black"), panel.grid.major = element_line(colour = "grey80"))
plot(perPenaltyCor_plot)







