#!/usr/bin/env Rscript
require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load("tidyverse")

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