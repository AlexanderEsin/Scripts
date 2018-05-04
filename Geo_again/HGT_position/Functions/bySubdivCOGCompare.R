#!/usr/bin/env Rscript
require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load("dplyr", "reshape2", "ggplot2", "ggdendro")

bySubdivCOGCompare	<- function(bySuvdiv_dataType_A, bySuvdiv_dataType_B, subDivision_list, clusterBy = NULL, subDivison_cols = NULL, cogMinCutOff = 50) {

	if (is.null(subDivison_cols)) stop("Required: subDivison_cols variable from the subDivisionKeyData object")

	# ------------------------------------------------------------------------------------- #

	combine_dataTypes			<- bind_rows(bySuvdiv_dataType_A, bySuvdiv_dataType_B)
	combine_dataTypes$subDiv	<- factor(combine_dataTypes$subDiv, levels = names(subDivision_list))

	# Remove all COG categories where there are fewer than 50 genes per COG in either dataType
	minCutOffCOGs	<- unlist(lapply(unique(combine_dataTypes$COGcat), function(COG) {
		sets	<- unique(combine_dataTypes$Set)
		subsets	<- unlist(lapply(sets, function(set) ifelse(sum(combine_dataTypes$numObsv[which(combine_dataTypes$Set == set & combine_dataTypes$COGcat == COG)]) >= cogMinCutOff, 1, 0)))
		result	<- ifelse(sum(subsets) != length(sets), NA, COG)
		return(result)
	}))
	minCutOffCOGs	<- minCutOffCOGs[!is.na(minCutOffCOGs)]
	minCutOffPrune	<- combine_dataTypes[which(combine_dataTypes$COGcat %in% minCutOffCOGs),]

	# ------------------------------------------------------------------------------------- #

	if (is.null(clusterBy)) {
		message("Clustering across all data")
		
		# Create a new column, combining COG names with the set (e.g. Ver) for clustering
		minCutOffPrune	<- bind_rows(lapply(unique(minCutOffPrune$Set), function(setName) {
			bySet_subdata	<- subset(minCutOffPrune, Set == setName)
			bySet_subdata$COGcat_wType	<- paste0(bySet_subdata$COGcat, "\n", setName)
			return(bySet_subdata)
		}))

		dataForCluster		<- minCutOffPrune
	} else if  (identical(clusterBy, "A")) {
		dataForCluster		<- minCutOffPrune[which(minCutOffPrune$Set == unique(bySuvdiv_dataType_A$Set)),]
	} else if (identical(clusterBy, "B")) {
		dataForCluster		<- minCutOffPrune[which(minCutOffPrune$Set == unique(bySuvdiv_dataType_B$Set)),]
	} else {
		stop("ClusterBy option should be either \'A\', \'B\', or \'NULL\'")
	}

	# ------------------------------------------------------------------------------------- #
	
	# Prepare data for clustering and calculate distance matrix
	if (is.null(clusterBy)) {
		clusterRecast_data				<- dcast(dataForCluster, COGcat_wType ~ subDiv, value.var = "numObsv")
		rownames(clusterRecast_data)	<- clusterRecast_data$COGcat_wType	
	} else {
		clusterRecast_data				<- dcast(dataForCluster, COGcat ~ subDiv, value.var = "numObsv")
		rownames(clusterRecast_data)	<- clusterRecast_data$COGcat
	}

	clusterRecastProp_data			<- sweep(clusterRecast_data[,-1], 1, rowSums(clusterRecast_data[,-1]), "/")
	clusterRecastProp_dist			<- dist(as.matrix(clusterRecastProp_data))

	# Produce dendrogram
	clusterCompartments_dendro		<- dendro_data(hclust(clusterRecastProp_dist, method = "ward.D2"))

	# Plot the clustering dendrogram
	perBranchCompartment_cluster	<- ggplot() +
		geom_segment(data = clusterCompartments_dendro$segments, aes(x = x, y = y, xend = xend, yend = yend), col = "#D9D9D9", size = 1) + 
		darkTheme +
		theme(
			plot.margin = unit(c(0, 1.4, 0, 1.4), "cm"),
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			axis.text = element_blank(),
			axis.title = element_blank(),
			axis.ticks = element_blank()
		)

	# ------------------------------------------------------------------------------------- #

	if (is.null(clusterBy)) {
		minCutOffPrune$COGcat_wType		<- factor(minCutOffPrune$COGcat_wType, levels = clusterCompartments_dendro$labels$label)

		compBySubDiv_barplot		<- ggplot(data = minCutOffPrune, aes(x = COGcat_wType, y = numObsv, fill = subDiv, label = as.character(numObsv))) +
			geom_bar(stat = "identity", position = "fill") +
			geom_text(position = position_fill(vjust = 0.5), color = "#333233") +
			# facet_wrap(~COGcat_wType, nrow = 1) +
			scale_fill_manual(values = subDivison_cols, guide = FALSE) +
			darkTheme +
			theme(
				plot.margin = unit(c(0, 2, 0.5, 2), "cm"),
				panel.grid.major.y = element_blank(),
				panel.grid.major.x = element_line(size = 0.8, color = "#D9D9D9"),
				panel.grid.minor.y = element_blank(),
				axis.text.y = element_blank(),
				axis.title.y = element_blank(),
				axis.title.x = element_blank(),
				axis.ticks = element_blank(),
				strip.background = element_rect(fill = "transparent", color = "#D9D9D9"),
				strip.text = element_text(color =  "#D9D9D9", size = 12)
			)

		# Return both the cluster dendrogram and barplot
		return(list(ClusterDendro = perBranchCompartment_cluster, ComparisonBarplot = compBySubDiv_barplot, perCOGstat = NA))

	} else {

		# Perform the ChisSQ tests - only if we are clustering by a type
		chisqByCOG_list	<- lapply(minCutOffCOGs, function(COG) {
			# subset data by COG
			subsetByCOG		<- minCutOffPrune[which(minCutOffPrune$COGcat == COG),]
			# recast data
			recastForChisq	<- dcast(subsetByCOG, subDiv ~ Set, value.var = "numObsv")
			# remove the category (subdiv) column leaving a two column matrix
			matrixForChisq	<- as.matrix(recastForChisq[,-1])
			# perform chisq
			chisq.pval		<- signif(chisq.test(x = matrixForChisq)$p.value, digits = 3)
			chisq.format	<- paste0("Chisq:\n", chisq.pval)
			return(data.frame(COGcat = COG, chisqPval = chisq.pval, chisqForm = chisq.format, stringsAsFactors = FALSE))
		})
		chisqByCOG_df		<- bind_rows(chisqByCOG_list)


		# Refactor the COG categories to match the clustering
		minCutOffPrune$COGcat		<- factor(minCutOffPrune$COGcat, levels = clusterCompartments_dendro$labels$label)
		chisqByCOG_df$COGcat		<- factor(chisqByCOG_df$COGcat, levels = clusterCompartments_dendro$labels$label)

		# Plot the comparison barplot
		compBySubDiv_barplot		<- ggplot(data = minCutOffPrune, aes(x = Set, y = numObsv, fill = subDiv, label = as.character(numObsv))) +
			geom_bar(stat = "identity", position = "fill") +
			geom_text(position = position_fill(vjust = 0.5), color = "#333233") +
			geom_label(data = chisqByCOG_df, aes(x = 1.5, y = 0.5, label = chisqForm), position = position_dodge(width = 0.5), size = 4, color = "#D9D9D9", fill = "#333233", inherit.aes = FALSE) +
			facet_wrap(~COGcat, nrow = 1) +
			scale_fill_manual(values = subDivison_cols, guide = FALSE) +
			darkTheme +
			theme(
				plot.margin = unit(c(0, 2, 0.5, 2), "cm"),
				panel.grid.major.y = element_blank(),
				panel.grid.major.x = element_line(size = 0.8, color = "#D9D9D9"),
				panel.grid.minor.y = element_blank(),
				axis.text.y = element_blank(),
				axis.title.y = element_blank(),
				axis.title.x = element_blank(),
				axis.ticks = element_blank(),
				strip.background = element_rect(fill = "transparent", color = "#D9D9D9"),
				strip.text = element_text(color =  "#D9D9D9", size = 12)
			)

		# Return both the cluster dendrogram and barplot
		return(list(ClusterDendro = perBranchCompartment_cluster, ComparisonBarplot = compBySubDiv_barplot, perCOGstat = chisqByCOG_df))

	}	
}

