#!/usr/bin/env Rscript

require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load("RSQLite", "dplyr", "ggplot2", "wesanderson", "parallel")

processOrthologPosition_2	<- function(COGname, orthData = NULL, dbConn = NULL, superGroupTaxid_list = NULL, dnaA_data = NULL, speciesSubgroups = NULL, outlierTaxid = NULL, distRange = c(0, 0.5)) {

	if (is.null(orthData)) stop("Provide an ortholog data set")
	if (is.null(dbConn)) stop("Provide a valid sqlite database connection")
	if (is.null(dnaA_data)) stop("Provide a data frame with taxids and their corresponding dnaA data")
	if (is.null(speciesSubgroups)) stop("Provide a data frame with species and corresponding subgroups (to account for species over-representation")
	if (is.null(superGroupTaxid_list)) stop("Provide superGroup breakdown") else superGroups	<- names(superGroupTaxid_list)

	message(paste0("COG == ", COGname))

	# First select dat based on the COG category
	onlyCOG_data		<- subset(orthData, COGcat == COGname)

	# Select groups that conform to the distRange cutoffs. For a group to be included, all genes within that group must lie within the distance range
	# I.e. if cut off = c(0, 0.25), exclude families where some genes are in cutoff and others are not
	inRangeGroups	<- lapply(unique(onlyCOG_data$orthGroup), function(group) {
		allGenes	<- subset(onlyCOG_data, orthGroup == group)
		allInRange	<- subset(onlyCOG_data, orthGroup == group & distToOri > min(distRange) & distToOri <= max(distRange))
		if(identical(nrow(allGenes), nrow(allInRange))) {
			return(group)
		} else {
			return(NA)
		}
	})
	inRangeGroups	<- sort(unlist(inRangeGroups[!is.na(inRangeGroups)]))

	# There might be no groups that satisfy the positional requirement in this COG
	if(length(inRangeGroups) == 0) {
		return(NA)
	}

	# Free up memory
	rm(list = c("orthData", "onlyCOG_data"))
	gc()

	byGroupOrthPosition	<- mclapply(inRangeGroups, function(group) {

		# All the data corresponding to the orthologs
		orthGroup_tbl	<- dbSendQuery(dbConn, 'SELECT * FROM t1 WHERE orthGroup = :orthGroup')
		dbBind(orthGroup_tbl, param = list(orthGroup = group))
		orthGroup_df	<- dbFetch(orthGroup_tbl)
		orthGroup_df	<- orthGroup_df[, !names(orthGroup_df) %in% c("sequence", "NuclSeq")]
		dbClearResult(orthGroup_tbl)

		# Select only those that have matching taxids in the clean DNAa_data
		pruneOrth_data	<- orthGroup_df[which(orthGroup_df$taxid %in% dnaA_data$taxid),]
		# Select only those that are on the main chromosome
		pruneOrth_data	<- pruneOrth_data[which(pruneOrth_data$plasmid == "F"),]

		# ------------------------------------------------------------------------------------- #
		# For each gene entry, combine with dnaA data (start, end, strand) and calculate relative gene position
		combineOri_df	<- left_join(pruneOrth_data, dnaA_clean_trim, by = "taxid")
		combineOri_df$relGeneStart	<- apply(combineOri_df, 1, genomeRelativePosition_format)
		combineOri_df$distToOri		<- ifelse(combineOri_df$relGeneStart > 0.5, 1 - combineOri_df$relGeneStart, combineOri_df$relGeneStart)

		# ------------------------------------------------------------------------------------- #
		# Separate AG genes from ortholog genes. We make no distinction about HGT or vertical here. Also remove the AG genomes that are rearranged
		combineOriNoAG_df	<- combineOri_df[which(combineOri_df$is_ag == 0),]
		combineOriAG_df		<- combineOri_df[which(combineOri_df$is_ag == 1 & !combineOri_df$taxid %in% outlierTaxid),]

		if (nrow(combineOriNoAG_df) == 0) {
			return(NA)
		}

		# ------------------------------------------------------------------------------------- #
		# Add subgrouping data
		combineOriSubgroup_df	<- left_join(combineOriNoAG_df, speciesSubgroups_df, by = "taxid")
		uniqueParentSpecies		<- unique(combineOriSubgroup_df$parentSpecies)
		numberOfGenes			<- nrow(combineOriSubgroup_df)

		bySubspeciesGroup_list	<- lapply(uniqueParentSpecies, summariseSubgroups, data = combineOriSubgroup_df, summaryStat = "distToOri")
		bySubspeciesGroup_df	<- bind_rows(bySubspeciesGroup_list)
		bySubspeciesGroup_df$Group	<- group
		
		# ------------------------------------------------------------------------------------- #
		bySuperGroup_list	<- lapply(superGroups, function(superGroup) {
			SGtaxid_list	<- superGroupTaxid_list[[superGroup]]
			subsetOrthPos	<- subset(bySubspeciesGroup_df, parentTaxid %in% SGtaxid_list)
			subsetOrthPos$SuperGroup	<- if (nrow(subsetOrthPos) > 0) superGroup else character(0)
			return(subsetOrthPos)
		})
		bySuperGroup	<- bind_rows(bySuperGroup_list)
		
		# ------------------------------------------------------------------------------------- #
		# Gene data for the HGT only for this group
		# combineOriAG_df$HGT	<- FALSE
		# combineOriAG_df$HGT[which(combineOriAG_df$protID %in% HGTOnlyData$protID)]	<- TRUE

		# ------------------------------------------------------------------------------------- #
		# # Make the plot
		# perGroupPlot	<- ggplot(data = bySubspeciesGroup_df, aes(x = distToOri, y = ..scaled..)) + 
		# 	scale_x_continuous(name = "Relative distance to Origin") +
		# 	scale_y_continuous(name = "Gene density of orthologs") +
		# 	geom_histogram(
		# 		data = subset(combineOriAG_df, HGT == FALSE),
		# 		aes(x = distToOri, y = ..count../nrow(combineOriAG_df)),
		# 		bins = 100, 
		# 		color = axisCol,
		# 		fill = wes_palette("Darjeeling2")[4],
		# 		alpha = 0.7,
		# 		inherit.aes = FALSE) +
		# 	geom_histogram(
		# 		data = subset(combineOriAG_df, HGT == TRUE),
		# 		aes(x = distToOri, y = ..count../nrow(combineOriAG_df)),
		# 		bins = 100, 
		# 		color = axisCol,
		# 		fill = wes_palette("Zissou1")[5],
		# 		alpha = 0.7,
		# 		inherit.aes = FALSE) +
		# 	stat_density(geom = "line", adjust = 1/10, color = axisCol, size = 0.8, alpha = 0.6) +
		# 	stat_density(data = subset(bySuperGroup, SuperGroup == "Bacillaceae"), geom = "line", adjust = 1/10, col = wes_palette("Zissou1")[3], size = 0.8, alpha = 0.9) +
		# 	ggtitle(paste0("COG = \'", COGname, "\' Gene family = ", group)) +
		# 	lightTheme

		# ------------------------------------------------------------------------------------- #
		# Clean memory
		rm(list = c("orthGroup_df", "combineOri_df", "combineOriNoAG_df", "combineOriSubgroup_df", "pruneOrth_data"))

		# return(list(allData = bySubspeciesGroup_df, bySuperGroupData = bySuperGroup, AG_HGT_data = combineOriAG_df[which(combineOriAG_df$HGT == TRUE),], numberOfGenes = numberOfGenes, numberOfSpecies = nrow(bySubspeciesGroup_df), plot = perGroupPlot))

		return(list(allData = bySubspeciesGroup_df, bySuperGroupData = bySuperGroup, numberOfGenes = numberOfGenes, numberOfSpecies = nrow(bySubspeciesGroup_df)))

	}, mc.cores = 10, mc.allow.recursive = FALSE)

	# Rename the list by group
	names(byGroupOrthPosition)	<- paste0(inRangeGroups)

	# ------------------------------------------------------------------------------------- #

	# Remove any groups with no orthologs
	byGroupOrthPosition		<- byGroupOrthPosition[!is.na(byGroupOrthPosition)]
	byGroupOrthPosition_df	<- bind_rows(lapply(byGroupOrthPosition, function(element) return(element$allData)))
	bySuperGroupCompare		<- bind_rows(lapply(byGroupOrthPosition, function(element) return(element$bySuperGroupData)))

	bySuperGroupCompare$SuperGroup		<- factor(bySuperGroupCompare$SuperGroup, levels = c(superGroups))
	bySuperGroupCompare$COGcat			<- COGname

	# bySuperGroup_withHGT			<- bind_rows(list(bySuperGroupCompare, HGTOnlyData[,c("distToOri", "SuperGroup")]))
	# bySuperGroup_withHGT$SuperGroup	<- factor(bySuperGroup_withHGT$SuperGroup, levels = c(superGroups, "HGT"))
	# bySuperGroup_withHGT$COGcat		<- COGname

	# # Set up colors
	# allCol	<- wes_palette("Chevalier1")[3]
	# grpCol	<- wes_palette("FantasticFox1")[1:length(superGroups)]
	# hgtCol	<- wes_palette("Darjeeling1")[1]

	# supergroupCompbyCOG_plot <- ggplot(data = byGroupOrthPosition_df, mapping = aes(distToOri)) +
	# 	scale_y_continuous(name = "Gene density") +
	# 	scale_x_continuous(name = "Relative distance from origin") +
	# 	stat_density(geom = "line", color = allCol, adjust = 1/10, position = "identity", size = 0.7) +
	# 	stat_density(data = bySuperGroup_withHGT, aes(distToOri, color = SuperGroup), geom = "line", adjust = 1/10, position = "identity", size = 0.7) +
	# 	scale_color_manual(values = c(grpCol, hgtCol), name = "Gene Group") +
	# 	ggtitle(paste0("COG = \'", COGname, "\'")) +
	# 	lightTheme + 
	# 	theme(
	# 		panel.grid.minor.y = element_blank(),
	# 		legend.background = element_rect(color = axisCol, fill = "white")
	# 	)

	# Clean memory
	rm(byGroupOrthPosition_df)
	gc()

	return(list(perGroupData = byGroupOrthPosition, supergroupCompareData = bySuperGroupCompare))
}
