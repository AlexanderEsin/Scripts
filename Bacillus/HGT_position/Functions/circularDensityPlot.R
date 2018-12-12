#!/usr/bin/env Rscript

require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load("circular", "wesanderson", "ggplot2", "polyclip")

circularDensityPlot	<- function(dataDensityA, dataDensityB = NULL, bgDensity, shrink = 1.4, titleCex = 1.2, titleName = "", tcl.offset = 1, uin = 2.1, bg = "white", axisCol = "black", enrichUpColor = NULL, enrichDownColor = NULL, densBLineColor = "black", axis.at = NULL, axis.labels = NULL) {

	# Check all data input is density.circular
	dataClass_l	<- lapply(list(dataDensityA, bgDensity), class)
	uniqueClass	<- unique(dataClass_l)
	if (length(uniqueClass) != 1 || uniqueClass != "density.circular") {
		stop("All data must be of class \'density.circular\'")
	}

	# Colour background
	if (!is.null(bg)) {
		par(bg = bg)
	}

	# Produce main plot and use the info for further plotting
	mainPlot <- circular:::plot.density.circular(
		x = bgDensity,
		points.plot = F,
		uin = uin,
		axes = F,
		xlab = NA,
		ylab = NA,
		control.circle = circle.control(type = "n"),
		nosort = T,
		lwd = 1.5,
		col = "transparent",
		ylim = c(-1.1, 1.1),
		plot.type = "circle",
		zero = (pi/2),
		rotation = "clock",
		shrink = shrink,
		main = ""
	)

	# Density lines
	bg_line	<- lines(bgDensity, lwd = 1.5, col = "transparent", plot.info = mainPlot, shrink = shrink)
	A_line	<- lines(dataDensityA, lwd = 0, col = "transparent", plot.info = mainPlot, shrink = shrink)

	# If colors are not provided, set up colors
	if (is.null(enrichUpColor))		enrichUpColor	<- wes_palette("Darjeeling1")[2]
	if (is.null(enrichDownColor))	enrichDownColor	<- wes_palette("Darjeeling1")[1]

	# Adjust the alpha of the colors for some transparency
	# enrichUpColor	<- alpha(enrichUpColor, 0.5)
	# enrichDownColor	<- alpha(enrichDownColor, 0.5)
	
	enrichDown	<- alpha(wes_palette("Darjeeling1")[1], 0.5)

	# Plot an invisible background line
	bg_line	<- lines(bgDensity, lwd = 1, col = "transparent", plot.info = mainPlot, shrink = shrink)

	# Colour polygons
	lapply(polyclip(A = list("x" = A_line$x, "y" = A_line$y), B = list("x" = bg_line$x, "y" = bg_line$y), op = "minus"), polygon, col = enrichUpColor, border = enrichUpColor)
	lapply(polyclip(A = list("x" = bg_line$x, "y" = bg_line$y), B = list("x" = A_line$x, "y" = A_line$y), op = "minus"), polygon, col = enrichDownColor, border = enrichDownColor)

	# If provided, plot the other data as a seperate line
	if (!is.null(dataDensityB) & identical(class(dataDensityB), "density.circular")) {
		B_line		<- lines(dataDensityB, lwd = 3.5, col = densBLineColor, plot.info = mainPlot, shrink = shrink)
	}

	# Replot the background line so it appears on top of the polygons
	bg_line	<- lines(bgDensity, lwd = 1, col = axisCol, plot.info = mainPlot, shrink = shrink)

	# Place the title in the middle
	text(0, 0.1, titleName, cex = titleCex, col = axisCol, font = 4)

	# Draw the axis - use custom functions to allow the the axis lines to be offset by a custom value
	par(cex.axis = titleCex)
	
	if (is.null(axis.at)) axis.labels <- c("Origin", "", "Terminus", "")
	else if (is.null(axis.at) & !is.null(axis.labels)) stop("If providing custom axis labels, provide custom axis coordinates with axis.at")
	else if (!is.null(axis.at) & is.null(axis.labels)) stop("If providing custom axis coordinates, provide custom axis labels with axis.labels")
	else if (!is.null(axis.at)) axis.at <- circular(axis.at, rotation = "clock", zero = (pi / 2), template = "none", units = "degrees")

	axis.circular_AE(at = axis.at, labels = axis.labels, rotation = "clock", zero = (pi / 2), template = "none", tcl = 0.3, tcl.text = 0.25, col = axisCol, tcl.offset = tcl.offset)

	# Record the plot as an object
	recordedPlot	<- recordPlot()
	return(recordedPlot)
}
