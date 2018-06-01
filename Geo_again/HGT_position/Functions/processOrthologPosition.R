#!/usr/bin/env Rscript

# Any libraries explicitly used in the script
require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load("RSQLite", "dplyr", "ggplot2", "wesanderson", "parallel")


processOrthologPosition	<- function(COGname, orthData = NULL, dbConn = NULL, superGroupTaxid_list = NULL, distRange = c(0, 0.5)) {

	if (is.null(orthData)) stop("Provide an ortholog data set")
	if (is.null(dbConn)) stop("Provide a valid sqlite database connection")
	if (is.null(superGroupTaxid_list)) stop("Provide superGroup breakdown")

	message(paste0("COG == ", COGname))

	# First select dat based on the COG category
	onlyCOG_data		<- subset(orthData, COGcat == COGname)

	# Select groups that conform to the distRange cutoffs. For a group to be included, all genes within that group must lie within the distance range
	# I.e. if cut off = c(0, 0.25), exclude families where some genes are in cutoff and others are not
	inRangeGroups	<- lapply(unique(onlyCOG_data$orthGroup), function(group) {
		allGenes	<- subset(onlyCOG_data, orthGroup == group)
		allInRange	<- subset(onlyCOG_data, orthGroup == group & dist2ori > min(distRange) & dist2ori <= max(distRange))
		if(identical(nrow(allGenes), nrow(allInRange))) {
			return(group)
		} else {
			return(NA)
		}
	})
	inRangeGroups	<- unlist(inRangeGroups[!is.na(inRangeGroups)])

	# There might be no groups that satisfy the positional requirement in this COG
	if(length(inRangeGroups) == 0) {
		return(NA)
	}

	# HGT only data to be combined later
	HGTOnlyData		<- onlyCOG_data[which(onlyCOG_data$orthGroup %in% inRangeGroups),]
	HGTOnlyData$SuperGroup	<- "HGT"

	byGroupOrthPosition	<- mclapply(inRangeGroups, function(group) {

		# All the data corresponding to the orthologs
		orthGroup_tbl	<- dbSendQuery(dbConn, 'SELECT * FROM t1 WHERE orthGroup = :orthGroup')
		dbBind(orthGroup_tbl, param = list(orthGroup = group))
		orthGroup_df	<- dbFetch(orthGroup_tbl)
		dbClearResult(orthGroup_tbl)

		# Select only those that have matching taxids in the clean DNAa_data
		pruneOrth_data	<- orthGroup_df[which(orthGroup_df$taxid %in% dnaA_clean_df$taxid),]
		# Select only those that are on the main chromosome
		pruneOrth_data	<- pruneOrth_data[which(pruneOrth_data$plasmid == "F"),]
		# Clean memory
		rm(orthGroup_df)

		# ------------------------------------------------------------------------------------- #
		# For each gene entry, combine with dnaA data (start, end, strand) and calculate relative gene position
		combineOri_df	<- left_join(pruneOrth_data, dnaA_clean_trim, by = "taxid")
		combineOri_df$relGeneStart	<- apply(combineOri_df, 1, genomeRelativePosition_format)
		combineOri_df$dist2ori		<- ifelse(combineOri_df$relGeneStart > 0.5, 1 - combineOri_df$relGeneStart, combineOri_df$relGeneStart)
		combineOri_df	<- subset(combineOri_df, select = -c(sequence, NuclSeq))

		# ------------------------------------------------------------------------------------- #
		# Separate AG genes from ortholog genes. We make no distinction about HGT or vertical here
		combineOriNoAG_df	<- combineOri_df[which(combineOri_df$is_ag == 0),]
		combineOriAG_df		<- combineOri_df[which(combineOri_df$is_ag == 1),]
		# Clean memory
		rm(combineOri_df)

		# ------------------------------------------------------------------------------------- #
		# Add subgrouping data
		combineOriSubgroup_df	<- left_join(combineOriNoAG_df, speciesSubgroups_df, by = "taxid")
		uniqueParentSpecies		<- unique(combineOriSubgroup_df$parentSpecies)

		# Clean memory
		rm(combineOriNoAG_df)

		bySubspeciesGroup_list	<- lapply(uniqueParentSpecies, summariseSubgroups, data = combineOriSubgroup_df, summaryStat = "dist2ori")
		bySubspeciesGroup_df	<- bind_rows(bySubspeciesGroup_list)
		
		# ------------------------------------------------------------------------------------- #
		# Gene data for the HGT only for this group
		HGT_entries		<- perTypeData$lHGT$`4`$allPosData$protID[which(perTypeData$lHGT$`4`$allPosData$orthGroup == group)]
		HGT_only_set	<- combineOriAG_df[which(combineOriAG_df$protID %in% HGT_entries),]


		# ------------------------------------------------------------------------------------- #
		# Make the plot
		plotTemp	<- ggplot(data = combineOriSubgroup_df, aes(dist2ori)) + 
			geom_density(adjust = 1/5) +
			geom_density(data = bySubspeciesGroup_df, adjust = 1/5, color = "red") +
			geom_histogram(
				data = combineOriAG_df,
				aes(x = dist2ori, y = ..ncount..),
				bins = 100, 
				fill = "red",
				inherit.aes = FALSE) +
			geom_histogram(
				data = HGT_only_set,
				aes(x = dist2ori, y = -..ncount..),
				bins = 100, 
				fill = "blue",
				inherit.aes = FALSE)

		gc()

		return(list(allData = bySubspeciesGroup_df, AG_HGT_data = HGT_only_set, numberOfGenes = nrow(combineOriSubgroup_df), numberOfSpecies = nrow(bySubspeciesGroup_df), plot = plotTemp))
	}, mc.cores = 10, mc.allow.recursive = FALSE)

	

	# Rename the list by group
	names(byGroupOrthPosition)	<- paste0(COGname, inRangeGroups)

	# ------------------------------------------------------------------------------------- #

	byGroupOrthPosition_df	<- bind_rows(lapply(byGroupOrthPosition, function(element) return(element$allData)))

	bySuperGroupCompare	<- bind_rows(lapply(superGroups, function(superGroup) {
		SGtaxid_list	<- superGroupTaxid_list[[superGroup]]
		subsetOrthPos	<- subset(byGroupOrthPosition_df, parentTaxid %in% SGtaxid_list)
		subsetOrthPos$SuperGroup	<- superGroup
		return(subsetOrthPos)
	}))

	bySuperGroupCompare	<- bind_rows(list(bySuperGroupCompare, HGTOnlyData))
	bySuperGroupCompare$SuperGroup	<- factor(bySuperGroupCompare$SuperGroup, levels = c(superGroups, "HGT"))
	bySuperGroupCompare$COGcat		<- COGname

	# Set up colors
	allCol	<- wes_palette("Chevalier1")[3]
	grpCol	<- wes_palette("FantasticFox1")[1:length(superGroups)]
	hgtCol	<- wes_palette("Darjeeling1")[1]

	supergroupCompbyCOG_plot <- ggplot(data = byGroupOrthPosition_df, mapping = aes(dist2ori)) +
		scale_y_continuous(name = "Gene density") +
		scale_x_continuous(name = "Relative distance from origin") +
		stat_density(geom = "line", color = allCol, adjust = 1/5, position = "identity") +
		stat_density(data = bySuperGroupCompare, aes(dist2ori, color = SuperGroup), geom = "line", adjust = 1/5, position = "identity") +
		scale_color_manual(values = c(grpCol, hgtCol), name = "Gene Group") +
		ggtitle(paste0("COG ", COGname)) +
		darkTheme + 
		theme(
			legend.background = element_rect(color = "#D9D9D9", fill = "#333233")
		)

	return(list(perGroupData = byGroupOrthPosition, supergroupCompareData = bySuperGroupCompare, supergroupComparePlot = supergroupCompbyCOG_plot))
}
