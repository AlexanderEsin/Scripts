#!/usr/bin/env Rscript

require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load("ggplot2", "GGally")

# GGpairs function to plot correlation subplots
cor_fun <- function(data, mapping, method = "pearson", ndp = 2, sz = 5, stars = TRUE, color = "red", ...) {

	data <- na.omit(data[,c(as.character(mapping$x)[2], as.character(mapping$y)[2])])

	x <- data[, as.character(mapping$x)[2]]
	y <- data[, as.character(mapping$y)[2]]

	corr	<- cor.test(x, y, method = method)
	est		<- corr$estimate

	if (stars) {
	  stars	<- c("***", "**", "*", "")[findInterval(corr$p.value, c(0, 0.001, 0.01, 0.05, 1))]
	  lbl	<- paste0(round(est, ndp), stars)
	} else {
	  lbl	<- round(est, ndp)
	}

	ggplot(data = data, mapping = mapping) +
	  annotate("text", x = mean(x), y = mean(y), label = lbl, size = sz, color = color, ...) +
	  theme(panel.grid = element_blank())
}


# GGpairs function to plot geom_smooth with points
smooth_lm_fun	<- function (data, mapping, smooth.colour = "black", method = "lm", ...) {
	p	<- ggplot(data = data, mapping)
	p	<- p + geom_point(...)
	if (!is.null(mapping$color) || !is.null(mapping$colour)) {
		p	<- p + geom_smooth(method = method)
	} else {
		p	<- p + geom_smooth(method = method, colour = smooth.colour, fill = alpha(smooth.colour, alpha = 0.05))
	}
	p
}