#!/usr/bin/env Rscript
require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load("tidyverse", "RColorBrewer")

# Load master variables and HGT position functions
invisible(sapply(HGTPos.all, source, .GlobalEnv))

# ------------------------------------------------------------------------------------- #
# Read in data

message("\nReading in data...", appendLF = FALSE)

# Per-type COG data
perTypeCOG_data		<- readRDS(file.path(positionData_path, "AG_perTypeCOGData.rds"))

# Zone list
zoneBoundaryList	<- readRDS(file.path(positionData_path, "AG_zoneBoundaries.rds"))

message("\rReading in data... done\n")

# ------------------------------------------------------------------------------------- #

# Output path for figures
COGEnrichment_path	<- file.path(figureOutput_path, "COG_enrichment")
if(!dir.exists(COGEnrichment_path)) dir.create(COGEnrichment_path)

# ------------------------------------------------------------------------------------- #

HGT_COGcount	<- perTypeCOG_data$lHGT$byZone %>%
	group_by(COGcat) %>%
	summarise(HGTtotalCount = sum(numObsv)) %>%
	mutate(HGTpropOfAll = HGTtotalCount / sum(HGTtotalCount))


Ver_COGcount	<- perTypeCOG_data$Ver$byZone %>%
	group_by(COGcat) %>%
	summarise(VerTotalCount = sum(numObsv)) %>%
	mutate(VerPropOfAll = VerTotalCount / sum(VerTotalCount))


HGTwithVer_COG	<- full_join(HGT_COGcount, Ver_COGcount, by = "COGcat") %>%
	mutate(HGT_enrichment = log(HGTpropOfAll / VerPropOfAll)) %>%
	arrange(desc(HGT_enrichment)) %>%
	mutate(COGcat = factor(COGcat, levels = COGcat)) %>%
	mutate(maxMin = case_when(
		HGT_enrichment > 0 ~ HGT_enrichment + 0.2,
		HGT_enrichment < 0 ~ HGT_enrichment - 0.2))

HGTwithVer_COG_N50	<- HGTwithVer_COG %>%
	subset(!HGTtotalCount < 50) %>%
	droplevels()

# ------------------------------------------------------------------------------------- #

HGT_COGenrichment_colorUse_barplot	<- ggplot(data = HGTwithVer_COG_N50, aes(x = COGcat, y = HGT_enrichment)) +
	scale_x_discrete(position = "top") +
	scale_y_continuous(
		name = "COG Enrichment in HGT genes",
		limits = c(-2.5, 2),
		breaks = seq(-2, 2, by = 1)) +
	geom_bar(aes(fill = HGT_enrichment), stat = "identity", color = axisCol) +
	scale_fill_gradientn(colours = c(dataTypeCols$Ver, "white", dataTypeCols$HGT), values = scales::rescale(c(min(HGTwithVer_COG_N50$HGT_enrichment), -0.5, 0, 0.5, max(HGTwithVer_COG_N50$HGT_enrichment)))) +
	# geom_label(aes(x = COGcat, y = (HGT_enrichment - 0) / 2,  label = HGTtotalCount), fill = "white") +
	geom_text(aes(x = COGcat, y = maxMin, label = COGcat), size = 6) +
	lightTheme + 
	theme(
		axis.ticks = element_blank(),
		panel.grid.minor.y = element_blank(),
		panel.grid.major.x = element_blank(),
		axis.line.x = element_blank(),
		axis.text.x = element_blank(),
		axis.title.x = element_blank()
	)



quartz(width = 20, height = 10)
print(HGT_COGenrichment_colorUse_barplot)
quartz.save(file = file.path(COGEnrichment_path, "HGT_v_Ver_COGEnrichment.pdf"), type = "pdf", dpi = 100)
invisible(dev.off())

#########

## Calculate across different HGT penalties ##
perTypeCOG_xPen_data <- lapply(penalty_list, function(hgtPenalty) {

	perTypeCOG_data <- lapply(dataTypes_withAge, function(dataType) {

		message(paste0("\tProcessing COG data for \'", dataType, "\'... "), appendLF = FALSE)

		byAge		<- FALSE
		subgroup	<- NA 

		if (identical(dataType, "All")) {
			penalty		<- NA
		} else if (identical(dataType, "Ver")) {
			penalty		<- verPenalty
		} else {
			penalty		<- hgtPenalty
			if (identical(dataType, "Old") || identical(dataType, "Recent")) {
				byAge		<- TRUE
				subgroup	<- ifelse(identical(dataType, "Old"), FALSE, TRUE)
			}
		}

		# For each COG, extract the relevant data - see the "processCOGData" function for exact outputs
		perCOGData_list			<- lapply(uniqueCOGs, processCOGData, data = perTypeData, dataType = dataType, penalty = penalty, byAge = byAge, subgroup = subgroup, zone_data = zoneBoundaryList$halfGenomeRange)
		names(perCOGData_list)	<- uniqueCOGs
		
		# Remove any COGs that don't have any observations
		perCOGData_list			<- perCOGData_list[!is.na(perCOGData_list)]

		# Produce the by-Compartment dataframe
		byCOGbyZone_df		<- bind_rows(lapply(perCOGData_list, function(COG) return(COG$perCOGbyZone_df)))
		byCOGbyZone_df$Set	<- dataType

		message(paste0("\r\tProcessing COG data for \'", dataType, penalty, "\'... done"))

		# Return all the raw datam as well as COGs per compartment, the circular summary of gene positions, and the circular weighted mean position plot
		return(list(perCOGdata = perCOGData_list, byZone = byCOGbyZone_df))
	})
	names(perTypeCOG_data)	<- dataTypes_withAge



	HGT_COGcount	<- perTypeCOG_data$Recent$byZone %>%
		group_by(COGcat) %>%
		summarise(HGTtotalCount = sum(numObsv)) %>%
		mutate(HGTpropOfAll = HGTtotalCount / sum(HGTtotalCount))

	print(sum(HGT_COGcount$HGTtotalCount))


	Ver_COGcount	<- perTypeCOG_data$Ver$byZone %>%
		group_by(COGcat) %>%
		summarise(VerTotalCount = sum(numObsv)) %>%
		mutate(VerPropOfAll = VerTotalCount / sum(VerTotalCount))

	print(sum(Ver_COGcount$VerTotalCount))


	HGTwithVer_COG	<- full_join(HGT_COGcount, Ver_COGcount, by = "COGcat") %>%
		mutate(HGT_enrichment = log(HGTpropOfAll / VerPropOfAll)) %>%
		arrange(desc(HGT_enrichment)) %>%
		mutate(COGcat = factor(COGcat, levels = COGcat)) %>%
		mutate(maxMin = case_when(
			HGT_enrichment > 0 ~ HGT_enrichment + 0.2,
			HGT_enrichment < 0 ~ HGT_enrichment - 0.2))

	HGTwithVer_COG_N50	<- HGTwithVer_COG %>%
		subset(!HGTtotalCount < 50) %>%
		droplevels()

	HGTwithVer_COG_N50 %<>% mutate(penalty = hgtPenalty)

	return(HGTwithVer_COG_N50)
})

x <- bind_rows(perTypeCOG_xPen_data)
x %<>% subset(!HGTtotalCount < 50) %>% mutate(COGcat = factor(COGcat, levels = levels(perTypeCOG_xPen_data[[1]]$COGcat))) %>% droplevels()

labelMinMax <- x %>% group_by(COGcat) %>% mutate(maxMinMax = case_when(maxMin < 0 ~ min(maxMin) - 0.1, maxMin > 0 ~ max(maxMin) + 0.1)) %>% filter(penalty == 3)

## Set up a colour ramp palette ##
col_ramp <- brewer.pal(9,"YlOrRd")
col_palette <- colorRampPalette(col_ramp[3:9])(4)

HGT_xPen_COGenrich_point	<- ggplot(data = x, aes(x = COGcat, y = HGT_enrichment, group = penalty, color = penalty)) +
	scale_x_discrete(position = "top") +
	scale_y_continuous(
		name = "COG Enrichment in HGT genes",
		limits = c(-2.5, 3),
		breaks = seq(-2, 3, by = 1)) +
	geom_point(stat = "identity", size = 3) +
	scale_color_manual(values = col_palette, name = "HGT Penalty") +
	geom_text(data = labelMinMax, aes(x = COGcat, y = maxMinMax, label = COGcat), size = 6, inherit.aes = FALSE) +
	lightTheme + 
	theme(
		axis.ticks = element_blank(),
		panel.grid.minor.y = element_blank(),
		panel.grid.major.x = element_blank(),
		axis.line.x = element_blank(),
		axis.text.x = element_blank(),
		axis.title.x = element_blank()
	)


coordinate_df <- data.frame(x1 = numeric(), x2 = numeric(), y1 = numeric(), y2 = numeric())

i = 1
for (COG in levels(x$COGcat)) {
	filter <- x %>% filter(COGcat == COG)

	## Calculate the max and min values for the boxes ##
	max_y <- max(filter$HGT_enrichment)
	min_y <- min(filter$HGT_enrichment)

	coordinate_df[nrow(coordinate_df)+1, ] <- c((i - 0.2), (i - 0.2), min_y, max_y)
	coordinate_df <- rbind(coordinate_df, c((i + 0.2), (i + 0.2), min_y, max_y))

	## Add horizontal lines ##
	coordinate_df <- rbind(coordinate_df, c((i - 0.2), (i + 0.2), min_y, min_y))
	coordinate_df <- rbind(coordinate_df, c((i - 0.2), (i + 0.2), max_y, max_y))

	i = i + 1
}

plot_relative_HGT_boxes <- HGT_xPen_COGenrich_point + geom_segment(data = coordinate_df, aes(x = x1, y = y1, xend = x2, yend = y2), inherit.aes = FALSE)

quartz(width = 12, height = 8)
print(plot_relative_HGT_boxes)
quartz.save(file = "/Users/aesin/Desktop/Thesis/CH3/Figs/Figure_3X_recentFunct_xfer.pdf", type = "pdf", dpi = 300)
invisible(dev.off())
