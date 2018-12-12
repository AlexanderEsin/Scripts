#!/usr/bin/env Rscript
require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load("dplyr", "ggplot2", "ggrepel", "wesanderson")

COGdistribCircular_plot	<- function(perCOGData_list, dataType, minCOGNum = 100, ymax = 0.7) {

	perCOGCircStart		<- lapply(perCOGData_list, function(COG) ifelse(length(COG$CircStart) >= minCOGNum, return(COG$CircStart), NA))
	perCOGCircTrim		<- perCOGCircStart[!is.na(perCOGCircStart)]


	# Apply the summary.circular_AE function to get, amongst other stats, the Mean and Rho values for each branch
	perCOGCircSummary_list	<- lapply(1:length(perCOGCircTrim), function(index) {
		# Get the branch, and summarise
		COGCircSum	<- summary.circular_AE(perCOGCircTrim[[index]])
		# Get the name of the branch
		COGName		<- names(perCOGCircTrim)[index]
		# Return data frame
		COGCircSum_df	<- cbind(COG = COGName, as.data.frame(t(COGCircSum)), stringsAsFactors = FALSE)
		return(COGCircSum_df)
	})
	# Bind the dataframes to make df with all stats per branch
	perCOGCircSummary_df		<- bind_rows(perCOGCircSummary_list)
	perCOGCircSummary_df$logN	<- log(perCOGCircSummary_df$n)


	conPal <- colorRampPalette(wes_palette("Zissou1"))

	# Produce the plot
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