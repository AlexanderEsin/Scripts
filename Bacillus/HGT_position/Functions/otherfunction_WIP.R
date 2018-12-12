





## // Subdivision functions // ##


## // Functional analysis // ##

## For a datatype and penalty (if not "All") extract COG list and total number
getPropCOG	<- function(dataType, Penalty = NA) {

	if (is.na(Penalty)) {
		COG			<- unlist(perTypeData[[dataType]]$allPosData$COGcat)
		colName		<- "Prop_All"
	} else {
		COG			<- unlist(perTypeData[[dataType]][[Penalty]]$allPosData$COGcat)
		colName		<- paste0(dataType, "_T", Penalty)
	}
	
	COG_df			<- as.data.frame(table(COG), stringsAsFactors = FALSE)
	totalCOGs		<- sum(COG_df$Freq)
	COG_df$Prop		<- COG_df$Freq / totalCOGs

	names(COG_df)[3]	<- colName
	return(list(propDF = COG_df[,-2], totalCOGs = totalCOGs))
}

## Calculate the enrichment of COGs in one Class (and penalty) over another
## For comparison against Vertical - use penalty 3 as this is the most stringent set
calcPropEnrich	<- function(data = "lHGT", background = "Ver", penalty = "3") {

	if (identical(data, background)) {
		stop("Can't compare to itself")
	} else if (identical(background, "Ver")) {
		bg_df	<- perTypeCOGprop$Ver$'3'$propCOGdf
	} else if (identical(background, "All")) {
		bg_df	<- perTypeCOGprop$All$propCOGdf
	}

	data_df		<- perTypeCOGprop[[data]][[penalty]]$propCOGdf

	join_df		<- inner_join(data_df, bg_df, by = "COG")
	join_df$Diff		<- log(join_df[,2] / join_df[,3])
	names(join_df)[4]	<- paste0("Diff_T", penalty)

	return(join_df[,c(1,4)])
}

## Calculate the per-COG ranges of enrichment values (across penalties)
## Produce coordinates to pass to geom_rect
calcEnrichRangesForPlot	<- function(df, molten_df, box_width = 0.2) {

	coordinate_df <- data.frame(x1 = numeric(), x2 = numeric(), y1 = numeric(), y2 = numeric())

	for (i in 1:length(df$COG)) {
		COG_needed	<- levels(molten_df$COG)[i]
		row_number	<- which(df$COG == COG_needed)

		## Do not count NAs in identifying min/max values ##
		data_line	<- df[row_number,-1]
		data_line_clean	<- data_line[,as.vector(!is.na(data_line))]

		## Calculate the max and min values for the boxes ##
		max_y	<- max(data_line_clean)
		min_y	<- min(data_line_clean)

		coordinate_df[nrow(coordinate_df)+1, ]	<- c((i - box_width), (i + box_width), min_y, max_y)
	}
	return(coordinate_df)
}

## // ggplot // ##
centerLegendTitle	<- function(plot) {
	## Requires gtools

	# extract legend
	g		<- ggplotGrob(plot)
	grobs 	<- g$grobs
	legend_index	<- which(sapply(grobs, function(x) x$name) == "guide-box")
	legend	<- grobs[[legend_index]]

	# extract guides table
	guides_index 	<- which(sapply(legend$grobs, function(x) x$name) == "layout")
	guides	<- legend$grobs[[guides_index]]

	# add extra column for spacing
	# guides$width[5] is the extra spacing from the end of the legend text
	# to the end of the legend title. If we instead distribute it 50:50 on
	# both sides, we get a centered legend
	guides	<- gtable_add_cols(guides, 0.5*guides$width[5], 1)
	guides$widths[6]	<- guides$widths[2]
	title_index	<- guides$layout$name == "title"
	guides$layout$l[title_index]	<- 2

	# reconstruct legend and write back
	legend$grobs[[guides_index]]	<- guides
	g$grobs[[legend_index]]			<- legend

	return(g)
}