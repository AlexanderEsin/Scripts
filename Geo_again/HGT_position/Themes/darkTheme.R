#!/usr/bin/env Rscript

require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load("ggplot2")

darkTheme	<- theme(
	plot.background = element_rect(fill = "#333233", color = NA),
	panel.background = element_rect(fill = "transparent", color = NA),
	panel.grid.major = element_line(size = 0.2, color = alpha("#D9D9D9", alpha = 0.5)),
	panel.grid.minor = element_line(size = 0.2, color = alpha("#D9D9D9", alpha = 0.5)),
	axis.ticks = element_line(size = 0.2, color = "#D9D9D9"),
	plot.title = element_text(size = 16, hjust = 0.5, color = "#D9D9D9"),
	legend.key = element_rect(fill = "transparent", color = NA),
	legend.background = element_rect(fill = "transparent", color = NA),
	legend.text = element_text(size = 14, color = "#D9D9D9"),
	legend.title = element_text(size = 14, color = "#D9D9D9"),
	legend.title.align = 0.5,
	axis.title = element_text(size = 14, color = "#D9D9D9"),
	axis.text = element_text(size = 14, colour = "#D9D9D9")
)