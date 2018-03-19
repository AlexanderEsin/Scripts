pointNumber	<- 50
wedgeSize	<- pi / 3
appearRange	<- pi / 2
pointFade	<- pi / 3


## Points axes (x = mean, y = rho)
pointsData_df	<- data.frame(Mean = runif(pointNumber, 0, 2*pi), Rho = runif(pointNumber, 0, 0.5))

# The values for the radar components - plot appearance etc is linked to the radarLineX variable
# Wedge is the trail, all Y values uniform
radar_df	<- data.frame(
	radarLineX = seq(0, 2*pi, length.out = 360),
	radarWedgeXstart = seq(0, 2*pi, length.out = 360),
	radarWedgeXend = seq(0, 2*pi, length.out = 360) - wedgeSize,
	radarAllYStart = rep(0, 360),
	radarAllYEnd = rep(0.5, 360),
	Frame = seq(1:360)
)

radarPointFade	<- runif(pointNumber, 0, pointFade)
radarPointAlpha	<- rescale(radarPointFade, to = c(0.2, 1))

## Account for negative wedge values (around the 0)
radar_df$radarWedgeXstart[which(radar_df$radarWedgeXend < 0)]	<- 2*pi + radar_df$radarWedgeXstart[which(radar_df$radarWedgeXend < 0)]
radar_df$radarWedgeXend[which(radar_df$radarWedgeXend < 0)]		<- 2*pi + radar_df$radarWedgeXend[which(radar_df$radarWedgeXend < 0)]

# Add the radarLine components to the points data (to calculate when points appear and disappear)
combined_list	<- lapply(1:nrow(radar_df), function(rowIndex) {
	row	<- radar_df[rowIndex,]
	pointsData_df$radarLineX		<- row$radarLineX
	pointsData_df$radarAllYStart	<- row$radarAllYStart
	pointsData_df$radarAllYEnd		<- row$radarAllYEnd
	pointsData_df$Frame				<- row$Frame

	pointsData_df$radarLineFade		<- unlist(lapply(1:nrow(pointsData_df), function(fadeIndex) {
		return(radarPointFade[fadeIndex])
	}))

	pointsData_df$Mean				<- ifelse((ifelse(pointsData_df$Mean + (appearRange + pointsData_df$radarLineFade) > 2*pi, pointsData_df$Mean - 2*pi, pointsData_df$Mean) < pointsData_df$radarLineX & ifelse(pointsData_df$Mean + (appearRange + pointsData_df$radarLineFade) > 2*pi, pointsData_df$Mean - 2*pi, pointsData_df$Mean) > pointsData_df$radarLineX - (appearRange + pointsData_df$radarLineFade)) | (pointsData_df$Mean < pointsData_df$radarLineX & pointsData_df$Mean > pointsData_df$radarLineX - (appearRange + pointsData_df$radarLineFade)), pointsData_df$Mean, NA)

	pointsData_df$Rho				<- ifelse((ifelse(pointsData_df$Mean + (appearRange + pointsData_df$radarLineFade) > 2*pi, pointsData_df$Mean - 2*pi, pointsData_df$Mean) < pointsData_df$radarLineX & ifelse(pointsData_df$Mean + (appearRange + pointsData_df$radarLineFade) > 2*pi, pointsData_df$Mean - 2*pi, pointsData_df$Mean) > pointsData_df$radarLineX - (appearRange + pointsData_df$radarLineFade)) | (pointsData_df$Mean < pointsData_df$radarLineX & pointsData_df$Mean > pointsData_df$radarLineX - (appearRange + pointsData_df$radarLineFade)), pointsData_df$Rho, NA)

	pointsData_df$radarPointAlphaVal		<- unlist(lapply(1:nrow(pointsData_df), function(fadeIndex) {
		radarPointAlphaVal	<- radarPointAlpha[fadeIndex]
		return(radarPointAlphaVal)
	}))
	return(pointsData_df)
})


for (frameIndex in 1:length(combined_list)) {
	if (frameIndex == 1) {
		combined_list[[frameIndex]]$radarPointCol	<- alpha(rep("green", pointNumber), combined_list[[frameIndex]]$radarPointAlphaVal)
		next
	}
	prevFrame			<- combined_list[[frameIndex - 1]]
	prevFramePlotted	<- which(!is.na(prevFrame$Mean))
	prevAlpha			<- prevFrame$radarPointAlphaVal[prevFramePlotted]
	newAlpha			<- prevAlpha * 0.95
	newAlpha[which(newAlpha < 0.05)]	<- 0.05
	combined_list[[frameIndex]]$radarPointAlphaVal[prevFramePlotted]	<- newAlpha
	combined_list[[frameIndex]]$radarPointCol							<- alpha(rep("green", pointNumber), combined_list[[frameIndex]]$radarPointAlphaVal)
	combined_list[[frameIndex]]$radarPointCol[prevFramePlotted]			<- alpha(rep("green", length(newAlpha)), newAlpha)
}
combined_df	<- bind_rows(combined_list)
combined_df	<- combined_df[!is.na(combined_df$Mean),]
nrow(combined_df[which(combined_df$Frame == 1),])



y <- lapply(seq(1:360), function(frame) {
	subset_df		<- combined_df[which(combined_df$Frame == frame),]

	plot	<- ggplot(subset_df, aes(x = Mean, y = Rho, color = radarPointCol)) +
		scale_x_continuous(labels = NULL, breaks = c(0, pi), limits = c(0, 2 * pi)) +
		scale_y_continuous(limits = c(0, 0.5)) +
		coord_polar() +
		# The background green over whole plot
		geom_rect(
			data = data.frame(xmin = 0, xmax = 2*pi, ymin = 0, ymax = 0.5),
			mapping = aes(
				xmin = xmin,
				xmax = xmax,
				ymin = ymin,
				ymax = ymax),
			inherit.aes = FALSE,
			fill = "darkgreen",
			color = NA,
			alpha = 0.4) +
		# The trailing wedge
		geom_rect(
			data = radar_df[frame,],
			mapping = aes(
				xmin = ifelse(radarWedgeXstart > 2*pi, radarWedgeXend, radarWedgeXend),
				xmax = ifelse(radarWedgeXstart > 2*pi, 2*pi, radarWedgeXstart),
				ymin = radarAllYStart,
				ymax = radarAllYEnd),
			inherit.aes = FALSE,
			fill = "green",
			color = NA,
			alpha = 0.2) +
		# If around 0, the other rect object (one for x>0 and one for x<0)
		geom_rect(
			data = radar_df[frame,],
			mapping = aes(
				xmin = ifelse(radarWedgeXstart > 2*pi, 0, 0),
				xmax = ifelse(radarWedgeXstart > 2*pi, radarWedgeXstart - 2*pi, 0),
				ymin = radarAllYStart,
				ymax = radarAllYEnd),
			inherit.aes = FALSE,
			fill = "green",
			color = NA,
			alpha = 0.2) +
		# Points plotted within a certain range following radar line 
		geom_point(
			data = subset_df,
			aes(
				x = ifelse(!is.na(Mean), Mean, NA),
				y = ifelse(!is.na(Mean), Rho, NA),
				colour = ifelse(!is.na(Mean), radarPointCol, NA)),
			size = 2.5,
			color = "green",
			inherit.aes = FALSE) +
		# Plot the radar leading line
		scale_color_manual(values = subset_df$radarPointCol, guide = FALSE) +
		geom_segment(data = subset_df, aes(x = radarLineX, xend = radarLineX, y = radarAllYStart, yend = radarAllYEnd), inherit.aes = FALSE, color = "green", size = 1.5) +
		# Theme
		theme(
			plot.background = element_rect(fill = "#333233", color = NA),
			panel.background = element_rect(fill = "transparent", color = NA),
			panel.grid = element_blank(),
			axis.ticks = element_blank(),
			plot.title = element_text(size = 16, hjust = 0.5, color = "#D9D9D9"),
			legend.key = element_rect(fill = "transparent", color = NA),
			legend.background = element_rect(fill = "transparent", color = NA),
			legend.text = element_text(size = 18, color = "#D9D9D9"),
			legend.title = element_text(size = 18, color = "#D9D9D9"),
			legend.title.align = 0.5,
			axis.title = element_text(size = 18, color = "#D9D9D9"),
			axis.text = element_text(size = 18, colour = "#D9D9D9")
		) +
		# The green axis lines have to be plotted independently
		geom_hline(yintercept = 0.1, color = "green", alpha = 0.8) +
		geom_hline(yintercept = 0.25, color = "green", alpha = 0.8) +
		geom_hline(yintercept = 0.5, color = "green", alpha = 0.8) +
		geom_vline(xintercept = 0, color = "green", alpha = 0.8) +
		geom_vline(xintercept = pi/2, color = "green", alpha = 0.8) +
		geom_vline(xintercept = pi, color = "green", alpha = 0.8) +
		geom_vline(xintercept = 1.5*pi, color = "green", alpha = 0.8)
		ggsave(plot = plot, filename = paste0("/Users/aesin/Desktop/Geo_again/ggRadar/Rplot", frame, ".png"), width = 6, height = 6, units = "in")
})


perBranchAvHGT_plot	<- ggplot(combined_df, aes(x = Mean, y = Rho, color = radarPointCol, frame = Frame)) +
	scale_x_continuous(labels = NULL, breaks = c(0, pi), limits = c(0, 2 * pi)) +
	scale_y_continuous(limits = c(0, 0.5)) +
	coord_polar() +
	# The background green over whole plot
	# geom_rect(
	# 	data = data.frame(xmin = 0, xmax = 2*pi, ymin = 0, ymax = 0.5),
	# 	mapping = aes(
	# 		xmin = xmin,
	# 		xmax = xmax,
	# 		ymin = ymin,
	# 		ymax = ymax),
	# 	inherit.aes = FALSE,
	# 	fill = "darkgreen",
	# 	color = NA,
	# 	alpha = 0.4) +
	# # The trailing wedge
	# geom_rect(
	# 	data = radar_df,
	# 	mapping = aes(
	# 		xmin = ifelse(radarWedgeXstart > 2*pi, radarWedgeXend, radarWedgeXend),
	# 		xmax = ifelse(radarWedgeXstart > 2*pi, 2*pi, radarWedgeXstart),
	# 		ymin = radarAllYStart,
	# 		ymax = radarAllYEnd,
	# 		frame = Frame),
	# 	inherit.aes = FALSE,
	# 	fill = "green",
	# 	color = NA,
	# 	alpha = 0.2) +
	# # If around 0, the other rect object (one for x>0 and one for x<0)
	# geom_rect(
	# 	data = radar_df,
	# 	mapping = aes(
	# 		xmin = ifelse(radarWedgeXstart > 2*pi, 0, 0),
	# 		xmax = ifelse(radarWedgeXstart > 2*pi, radarWedgeXstart - 2*pi, 0),
	# 		ymin = radarAllYStart,
	# 		ymax = radarAllYEnd,
	# 		frame = Frame),
	# 	inherit.aes = FALSE,
	# 	fill = "green",
	# 	color = NA,
	# 	alpha = 0.2) +
	# Points plotted within a certain range following radar line 
	geom_point(
		data = combined_df,
		aes(
			# x = ifelse((ifelse(Mean + (appearRange + radarLineFade) > 2*pi, Mean - 2*pi, Mean) < radarLineX & ifelse(Mean + (appearRange + radarLineFade) > 2*pi, Mean - 2*pi, Mean) > radarLineX - (appearRange + radarLineFade)) | (Mean < radarLineX & Mean > radarLineX - (appearRange + radarLineFade)), Mean, NA),
			# y = ifelse((ifelse(Mean + (appearRange + radarLineFade) > 2*pi, Mean - 2*pi, Mean) < radarLineX & ifelse(Mean + (appearRange + radarLineFade) > 2*pi, Mean - 2*pi, Mean) > radarLineX - (appearRange + radarLineFade)) | (Mean < radarLineX & Mean > radarLineX - (appearRange + radarLineFade)), Rho, NA),
			# colour = ifelse((ifelse(Mean + (appearRange + radarLineFade) > 2*pi, Mean - 2*pi, Mean) < radarLineX & ifelse(Mean + (appearRange + radarLineFade) > 2*pi, Mean - 2*pi, Mean) > radarLineX - (appearRange + radarLineFade)) | (Mean < radarLineX & Mean > radarLineX - (appearRange + radarLineFade)), radarPointCol, NA),
			x = ifelse(!is.na(Mean), Mean, NA),
			y = ifelse(!is.na(Mean), Rho, NA),
			colour = ifelse(!is.na(Mean), radarPointCol, NA),
			frame = Frame),
		size = 2.5,
		# color = "green",
		inherit.aes = FALSE) +
	# Plot the radar leading line
	# guides(color = FALSE) +
	# scale_color_manual() +
	geom_segment(data = combined_df, aes(x = radarLineX, xend = radarLineX, y = radarAllYStart, yend = radarAllYEnd, frame = Frame), inherit.aes = FALSE, color = "green", size = 1.5) +
	# Theme
	theme(
		plot.background = element_rect(fill = "#333233", color = NA),
		panel.background = element_rect(fill = "transparent", color = NA),
		panel.grid = element_blank(),
		axis.ticks = element_blank(),
		plot.title = element_text(size = 16, hjust = 0.5, color = "#D9D9D9"),
		legend.key = element_rect(fill = "transparent", color = NA),
		legend.background = element_rect(fill = "transparent", color = NA),
		legend.text = element_text(size = 18, color = "#D9D9D9"),
		legend.title = element_text(size = 18, color = "#D9D9D9"),
		legend.title.align = 0.5,
		axis.title = element_text(size = 18, color = "#D9D9D9"),
		axis.text = element_text(size = 18, colour = "#D9D9D9")
	) +
	# The green axis lines have to be plotted independently
	geom_hline(yintercept = 0.1, color = "green", alpha = 0.8) +
	geom_hline(yintercept = 0.25, color = "green", alpha = 0.8) +
	geom_hline(yintercept = 0.5, color = "green", alpha = 0.8) +
	geom_vline(xintercept = 0, color = "green", alpha = 0.8) +
	geom_vline(xintercept = pi/2, color = "green", alpha = 0.8) +
	geom_vline(xintercept = pi, color = "green", alpha = 0.8) +
	geom_vline(xintercept = 1.5*pi, color = "green", alpha = 0.8)

gganimate(perBranchAvHGT_plot, ani.width = 1000, ani.height = 1000, ani.res = 1200, interval = 0.025, title_frame = FALSE, filename = "/Users/aesin/Desktop/RadarRand.mp4")