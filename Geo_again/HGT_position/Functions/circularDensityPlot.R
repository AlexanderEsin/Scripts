#!/usr/bin/env Rscript

require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load("circular", "wesanderson", "ggplot2", "polyclip")

circularDensityPlot	<- function(dataDensityA, dataDensityB = NULL, bgDensity, shrink = 1.4, titleCex = 1.2, titleName = "", tcl.offset = 1, uin = 2.1, bg = "#333233") {

	# Check all data input is density.circular
	dataClass_l	<- lapply(list(dataDensityA, bgDensity), class)
	uniqueClass	<- unique(dataClass_l)
	if (length(uniqueClass) != 1 || uniqueClass != "density.circular") {
		stop("All data must be of class \'density.circular\'")
	}

	# Colour background
	if (!is.na(bg)) {
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

	# Set up colours
	enrichUp	<- alpha(wes_palette("Darjeeling1")[2], 0.5)
	enrichDown	<- alpha(wes_palette("Darjeeling1")[1], 0.5)

	# Plot an invisible background line
	bg_line	<- lines(bgDensity, lwd = 1, col = "transparent", plot.info = mainPlot, shrink = shrink)

	# Colour polygons
	lapply(polyclip(A = list("x" = A_line$x, "y" = A_line$y), B = list("x" = bg_line$x, "y" = bg_line$y), op = "minus"), polygon, col = enrichUp, border = enrichUp)
	lapply(polyclip(A = list("x" = bg_line$x, "y" = bg_line$y), B = list("x" = A_line$x, "y" = A_line$y), op = "minus"), polygon, col = enrichDown, border = enrichDown)

	# If provided, plot the other data as a seperate line
	if (!is.null(dataDensityB) & identical(class(dataDensityB), "density.circular")) {
		Bline_col	<- alpha(wes_palette("Darjeeling1")[3], 0.8)
		B_line		<- lines(dataDensityB, lwd = 1.5, col = Bline_col, plot.info = mainPlot, shrink = shrink)
	}

	# Replot the background line so it appears on top of the polygons
	bg_line	<- lines(bgDensity, lwd = 1, col = "#D9D9D9", plot.info = mainPlot, shrink = shrink)

	# Place the title in the middle
	text(0, 0.2, titleName, cex = titleCex, col = "#D9D9D9", font = 4)

	# Draw the axis - use custom functions to allow the the axis lines to be offset by a custom value
	par(cex.axis = titleCex)
	axis.circular_AE(at = NULL, labels = c("Origin", "", "Terminus", ""), rotation = "clock", zero = (pi / 2), template = "none", tcl = 0.3, tcl.text = 0.25, col = "#D9D9D9", tcl.offset = tcl.offset)

	# Record the plot as an object
	recordedPlot	<- recordPlot()
	return(recordedPlot)
}
