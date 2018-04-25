summariseSubgroups	<- function(df, subgroupData, summaryStat, clusterCon = NULL) {
	
	# Check there is a cluster
	if (is.null(clusterCon)) {
		stop("Need to specify a cluster connection for this function")
	}

	uniqueParentTaxids	<- unique(subgroupData$internalParent)

	byParentTaxid	<- parLapply(clusterCon, uniqueParentTaxids, function(parentTaxid) {
		# Get all the data rows that correspond to children of this parent
		childTaxids		<- subgroupData$taxid[which(subgroupData$internalParent == parentTaxid)]
		perParentData	<- df[which(df$taxid %in% childTaxids),]

		# Some taxids will not be present in subsetted data
		if(nrow(perParentData) == 0) {
			return(NA)
		}

		# The parent taxid might not be present in the data, so we take the name from the subgrouping data
		parentName		<- subgroupData$parentName[which(subgroupData$taxid == parentTaxid)]

		# Summary statistics on the column desired
		byParent_mean	<- mean(perParentData[[summaryStat]])
		byParent_sd		<- ifelse(nrow(perParentData) > 1, sd(perParentData[[summaryStat]]), NA)

		# Make and rename output df
		out_df	<- data.frame(parentTaxid = parentTaxid, parentBinomial = parentName, childNumber = nrow(perParentData), summaryMean = byParent_mean, summarySD = byParent_sd, stringsAsFactors = FALSE)
		names(out_df)[which(names(out_df) == "summaryMean")]	<- summaryStat

		return(out_df)
	})
	# Remove empty element and bind to df. Then order by number of children
	byParentTaxid_df	<- bind_rows(byParentTaxid[!is.na(byParentTaxid)])
	byParentTaxid_df	<- byParentTaxid_df[order(byParentTaxid_df$childNumber, decreasing = TRUE),]
	return(byParentTaxid_df)
}