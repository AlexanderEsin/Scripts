#!/usr/bin/env Rscript

# Load master variables and HGT position functions
invisible(sapply(HGTPos.all, source, .GlobalEnv))

# Any libraries explicitly used in the script
require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load("dplyr", "reshape2", "grid", "ggplot2", "ggrepel", "ggdendro", "ggpubr", "wesanderson", "RColorBrewer")

# ------------------------------------------------------------------------------------- #
message("\nReading in data...", appendLF = FALSE)

# General position data
perTypeData			<- readRDS(file.path(positionData_path, "AG_perTypeData.rds"))

# Zone boundaries
zoneBoundaryList	<- readRDS(file.path(positionData_path, "AG_zoneBoundaries.rds"))

message("\rReading in data... done\n")


HGTdata	<- perTypeData$lHGT$'4'$allPosData %>% select(-c(CircStart, CircEnd)) %>% mutate(type = "HGT")
Verdata	<- perTypeData$Ver$'3'$allPosData %>% select(-c(CircStart, CircEnd)) %>% mutate(type = "Ver")
Alldata	<- perTypeData$All$allPosData %>%  select(-c(CircStart, CircEnd)) %>% mutate(type = "All")

# Add the age - seperated data together
HGT_byAge	<- HGTdata %>%
	mutate(type = case_when(
		Subgroup == TRUE ~ "Recent",
		TRUE ~ "Old")) %>%
	bind_rows(HGTdata, .)

# Combine vertical and HGT genes
combined_HGTVer	<- HGT_byAge %>% bind_rows(Verdata)

# The "gray" zone is the "all" data excluding all vertical and HGT genes
combined_all	<- Alldata %>% 
	subset(!protID %in% combined_HGTVer$protID) %>%
	mutate(type = "Other") %>%
	bind_rows(combined_HGTVer, .) %>%
	mutate(type = factor(type, levels = c("HGT", "Other", "Ver", "Recent", "Old"))) %>%
	mutate(plotDiv = case_when(
		type == "Recent" | type == "Old" ~ "byAge",
		TRUE ~ "byType")) %>% 
	mutate(plotDiv = factor(plotDiv, levels = c("byType", "byAge"))) %>%
	group_by(type, plotDiv) %>%
	summarise(numberOfGenes = n())



countGenesByType_barplot	<- ggplot(data = combined_all, aes(x = plotDiv, y = numberOfGenes, label = as.character(numberOfGenes))) +
	geom_bar(aes(fill = type), stat = "identity", position = "fill") +
	geom_label(aes(group = type), position = position_fill(vjust = 0.5), fill = "white", show.legend = FALSE) +
	scale_fill_manual(values = c(dataTypeCols$HGT, dataTypeCols$Other, dataTypeCols$Ver, dataTypeCols$Recent, dataTypeCols$Old)) +
	lightTheme


quartz(width = 12, height = 8)
print(countGenesByType_barplot)
quartz.save(file = file.path(figureOutput_path, "geneNumberStats.pdf"), type = "pdf", dpi = 300)
invisible(dev.off())