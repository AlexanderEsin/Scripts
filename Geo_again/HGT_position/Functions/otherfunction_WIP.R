





## // Subdivision functions // ##



plotCOGCircularDistribution	<- function(perCOGData_list, dataType, minCOGNum = 100, ymax = 0.7) {

	perCOGCircStart		<- lapply(perCOGData_list, function(COG) ifelse(length(COG$CircStart) >= minCOGNum, return(COG$CircStart), NA))
	perCOGCircTrim		<- perCOGCircStart[!is.na(perCOGCircStart)]


	## Apply the summary.circular_AE function to get, amongst other stats, the Mean and Rho values for each branch
	perCOGCircSummary_list	<- lapply(1:length(perCOGCircTrim), function(index) {
		# Get the branch, and summarise
		COGCircSum	<- summary.circular_AE(perCOGCircTrim[[index]])
		# Get the name of the branch
		COGName		<- names(perCOGCircTrim)[index]
		# Return data frame
		COGCircSum_df	<- cbind(COG = COGName, as.data.frame(t(COGCircSum)), stringsAsFactors = FALSE)
		return(COGCircSum_df)
	})
	## Bind the dataframes to make df with all stats per branch
	perCOGCircSummary_df		<- bind_rows(perCOGCircSummary_list)
	perCOGCircSummary_df$logN	<- log(perCOGCircSummary_df$n)


	conPal <- colorRampPalette(wes_palette("Zissou1"))
	perCOGAvPos_plot	<- ggplot(perCOGCircSummary_df, aes(x = Mean, y = Rho, color = logN, label = COG)) +
		scale_x_continuous(labels = c("Origin", "Terminus"), breaks = c(0, pi), limits = c(0, 2 * pi)) +
		scale_y_continuous(limits = c(0, ymax)) +
		ggtitle(paste0("Average COG position for ", dataType, " genes")) +
		coord_polar() +
		geom_point(size = 1.5) +
		scale_colour_gradientn(colors = conPal(20)) +
		geom_text_repel(aes(label = COG), point.padding = 0.2, min.segment.length = 0.2) +
		darkTheme +
		theme(
			axis.title.x = element_blank(),
			legend.justification = c(0, 0),
			legend.position = c(0, 0),
			panel.grid.major = element_line(size = 0.2, color = alpha("#D9D9D9", alpha = 0.3)),
			panel.grid.minor = element_line(size = 0.2, color = alpha("#D9D9D9", alpha = 0.3)),
		)

	return(list(circSum = perCOGCircSummary_df, plot = perCOGAvPos_plot))
}

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

plotBranchbyTrait_AE <- function (tree, x, mode = c("edges", "tips", "nodes"), palette = "rainbow", legend = TRUE, xlims = NULL, ...) {
	
	mode <- mode[1]
	if (!inherits(tree, "phylo")) 
		stop("tree should be an object of class \"phylo\".")
	if (mode == "tips") {
		x <- c(x[tree$tip.label], fastAnc(tree, x))
		names(x)[1:length(tree$tip.label)] <- 1:length(tree$tip.label)
		XX <- matrix(x[tree$edge], nrow(tree$edge), 2)
		x <- rowMeans(XX)
	}
	else if (mode == "nodes") {
		XX <- matrix(x[tree$edge], nrow(tree$edge), 2)
		x <- rowMeans(XX)
	}
	if (hasArg(tol)) 
		tol <- list(...)$tol
	else tol <- 1e-06
	if (hasArg(prompt)) 
		prompt <- list(...)$prompt
	else prompt <- FALSE
	if (hasArg(type)) 
		type <- list(...)$type
	else type <- "phylogram"
	if (hasArg(show.tip.label)) 
		show.tip.label <- list(...)$show.tip.label
	else show.tip.label <- TRUE
	if (hasArg(show.node.label)) 
		show.node.label <- list(...)$show.node.label
	else show.node.label <- FALSE
	if (hasArg(edge.width)) 
		edge.width <- list(...)$edge.width
	else edge.width <- 4
	if (hasArg(edge.lty)) 
		edge.lty <- list(...)$edge.lty
	else edge.lty <- 1
	if (hasArg(font)) 
		font <- list(...)$font
	else font <- 3
	if (hasArg(cex)) 
		cex <- list(...)$cex
	else cex <- par("cex")
	if (hasArg(adj)) 
		adj <- list(...)$adj
	else adj <- NULL
	if (hasArg(srt)) 
		srt <- list(...)$srt
	else srt <- 0
	if (hasArg(no.margin)) 
		no.margin <- list(...)$no.margin
	else no.margin <- TRUE
	if (hasArg(root.edge)) 
		root.edge <- list(...)$root.edge
	else root.edge <- FALSE
	if (hasArg(label.offset)) 
		label.offset <- list(...)$label.offset
	else label.offset <- 0.01 * max(nodeHeights(tree))
	if (hasArg(underscore)) 
		underscore <- list(...)$underscore
	else underscore <- FALSE
	if (hasArg(x.lim)) 
		x.lim <- list(...)$x.lim
	else x.lim <- NULL
	if (hasArg(y.lim)) 
		y.lim <- list(...)$y.lim
	else y.lim <- if (legend && !prompt && type %in% c("phylogram", 
		"cladogram")) 
		c(1 - 0.06 * length(tree$tip.label), length(tree$tip.label))
	else NULL
	if (hasArg(direction)) 
		direction <- list(...)$direction
	else direction <- "rightwards"
	if (hasArg(lab4ut)) 
		lab4ut <- list(...)$lab4ut
	else lab4ut <- NULL
	if (hasArg(tip.color)) 
		tip.color <- list(...)$tip.color
	else tip.color <- "black"
	if (hasArg(plot)) 
		plot <- list(...)$plot
	else plot <- TRUE
	if (hasArg(rotate.tree)) 
		rotate.tree <- list(...)$rotate.tree
	else rotate.tree <- 0
	if (hasArg(open.angle)) 
		open.angle <- list(...)$open.angle
	else open.angle <- 0

	if (is.function(palette)) 
		cols <- palette(n = 1000)
	else {
		if (palette == "heat.colors") 
			cols <- heat.colors(n = 1000)
		if (palette == "gray") 
			cols <- gray(1000:1/1000)
		if (palette == "rainbow") 
			cols <- rainbow(1000, start = 0.7, end = 0)
		else
			cols <- palette
	}
	if (is.null(xlims)) 
		xlims <- range(x) + c(-tol, tol)
	breaks <- 0:1000/1000 * (xlims[2] - xlims[1]) + xlims[1]
	whichColor <- function(p, cols, breaks) {
		i <- 1
		while (p >= breaks[i] && p > breaks[i + 1]) i <- i + 
			1
		cols[i]
	}
	colors <- sapply(x, whichColor, cols = cols, breaks = breaks)
	par(lend = 2)
	xx <- plot.phylo(tree, type = type, show.tip.label = show.tip.label, 
		show.node.label = show.node.label, edge.color = colors, 
		edge.width = edge.width, edge.lty = edge.lty, font = font, 
		cex = cex, adj = adj, srt = srt, no.margin = no.margin, 
		root.edge = root.edge, label.offset = label.offset, underscore = underscore, 
		x.lim = x.lim, y.lim = y.lim, direction = direction, 
		lab4ut = lab4ut, tip.color = tip.color, plot = plot, 
		rotate.tree = rotate.tree, open.angle = open.angle, lend = 2, 
		new = FALSE)
	
	if (legend == TRUE && is.logical(legend)) 
		legend <- round(0.3 * max(nodeHeights(tree)), 2)
	if (legend) {
		if (hasArg(title)) 
			title <- list(...)$title
		else title <- "trait value"
		if (hasArg(digits)) 
			digits <- list(...)$digits
		else digits <- 1
		if (prompt) 
			add.color.bar(legend, cols, title, xlims, digits, 
				prompt = TRUE)
		else add.color.bar(legend, cols, title, xlims, digits, 
			prompt = FALSE, x = par()$usr[1] + 0.05 * (par()$usr[2] - 
				par()$usr[1]), y = par()$usr[3] + 0.05 * (par()$usr[4] - 
				par()$usr[3]))
	}
	invisible(xx)
}
