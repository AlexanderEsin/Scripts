#!/usr/bin/env Rscript

summariseSubgroups	<- function(parentTaxid, data = NULL, summaryStat = "dist2ori") {

	if (is.null(data)) stop("Need to provide dataframe")

	# Subset by the parent species taxid
	byParent_data	<- data[which(data$parentSpecies == parentTaxid),]

	# Binomial name of the parent species
	parentName		<- unique(byParent_data$parentName)

	# Summary statistics on the column desired
	byParent_mean	<- mean(byParent_data[[summaryStat]])
	byParent_sd		<- ifelse(nrow(byParent_data) > 1, sd(byParent_data[[summaryStat]]), NA)

	# Make and rename output df
	out_df	<- data.frame(
		parentTaxid = parentTaxid,
		parentBinomial = parentName,
		childNumber = nrow(byParent_data),
		summaryMean = byParent_mean,
		summarySD = byParent_sd,
		stringsAsFactors = FALSE)

	# Rename the "mean" column to the name of the column which we averaged
	names(out_df)[which(names(out_df) == "summaryMean")]	<- summaryStat

	return(out_df)
}