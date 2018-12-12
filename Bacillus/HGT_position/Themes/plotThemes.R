#!/usr/bin/env Rscript

require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load("ggplot2", "wesanderson")

axisCol		<- wes_palette("Moonrise2")[4]
textCol		<- wes_palette("Moonrise2")[4]
lightTheme	<- theme(
	plot.background = element_rect(fill = "white", color = NA),
	panel.background = element_rect(fill = "transparent", color = NA),
	panel.grid.major = element_line(size = 0.4, color = alpha(axisCol, alpha = 0.5)),
	panel.grid.minor = element_line(size = 0.2, color = alpha(axisCol, alpha = 0.5)),
	axis.ticks = element_line(size = 0.2, color = axisCol),
	plot.title = element_text(size = 16, hjust = 0.5, color = textCol),
	legend.key = element_rect(fill = "transparent", color = NA),
	legend.background = element_rect(fill = "transparent", color = NA),
	legend.text = element_text(size = 14, color = textCol),
	legend.title = element_text(size = 14, color = textCol),
	legend.title.align = 0.5,
	axis.title = element_text(size = 14, color = textCol),
	axis.text = element_text(size = 14, colour = textCol)
)


darkAxisCol	<- "#D9D9D9"
darkTextCol	<- "#D9D9D9"

darkTheme	<- theme(
	plot.background = element_rect(fill = "#333233", color = NA),
	panel.background = element_rect(fill = "transparent", color = NA),
	panel.grid.major = element_line(size = 0.2, color = alpha(darkAxisCol, alpha = 0.5)),
	panel.grid.minor = element_line(size = 0.2, color = alpha(darkAxisCol, alpha = 0.5)),
	axis.ticks = element_line(size = 0.2, color = darkAxisCol),
	plot.title = element_text(size = 16, hjust = 0.5, color = darkTextCol),
	legend.key = element_rect(fill = "transparent", color = NA),
	legend.background = element_rect(fill = "transparent", color = NA),
	legend.text = element_text(size = 14, color = darkTextCol),
	legend.title = element_text(size = 14, color = darkTextCol),
	legend.title.align = 0.5,
	axis.title = element_text(size = 14, color = darkTextCol),
	axis.text = element_text(size = 14, colour = darkTextCol)
)