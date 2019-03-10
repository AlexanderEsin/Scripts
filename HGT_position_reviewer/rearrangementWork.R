invisible(sapply(HGTPos.geo, source, .GlobalEnv))

require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load("tidyverse", "reshape2", "grid", "ggplot2", "ggrepel", "ggdendro", "ggpubr", "wesanderson", "RColorBrewer", "fs")

out_path <- "/Users/aesin/Documents/Scripts/HGT_position_reviewer/Figures/rearrangement"
dir_create(out_path)

# ------------------------------------------------------------------------------------- #
message("\nReading in data...", appendLF = FALSE)

# General position data
GPA_perTypeData			<- readRDS(file.path(positionData_path, "AG_perTypeData_25genomes.rds"))

# Zones
GPA_zoneBoundaryList	<- readRDS(file.path(positionData_path, "AG_zoneBoundaries.rds"))

# ------------------------------------------------------------------------------------- #

Y412M_taxids <- c(544556, 550542)

labelledY412M_data <- GPA_perTypeData$All$allPosData %>%
	filter(taxid %in% Y412M_taxids) %>%
	select(-c(CircStart, CircEnd)) %>%
	mutate(Type = case_when(
		protID %in% GPA_perTypeData$lHGT$`4`$allPosData$protID ~ "HGT",
		protID %in% GPA_perTypeData$Ver$`3`$allPosData$protID ~ "Ver",
		TRUE ~ "Other"))

withExtend <- labelledY412M_data %>%
	filter(relGeneStart < 0.2 | relGeneStart > 0.8) %>%
	mutate(relGeneStart = case_when(
		relGeneStart < 0.2 ~ relGeneStart + 1,
		TRUE ~ relGeneStart - 1)
	) %>%
	bind_rows(labelledY412M_data)


# ------------------------------------------------------------------------------------- #

mc_52_max <- labelledY412M_data %>% filter(taxid == 550542 & relGeneEnd == max(relGeneEnd)) %>% pull(geneEnd)
mauveBounds_mc52 <- tibble(block = rep(seq(1:5), each  = 2),
					  orient = rep(c("start", "end"), 5),
					  coord = c(0, 170000, 170001, 2190000, 2190001, 2567000, 2567001, 3060000, 3060001, 3628000))

mc_61_max <- labelledY412M_data %>% filter(taxid == 544556) %>% filter(relGeneEnd == max(relGeneEnd)) %>% pull(geneEnd)
mauveBounds_mc61 <- tibble(block = rep(seq(1:5), each  = 2),
						   orient = rep(c("start", "end"), 5),
						   coord = c(0, 170000, 1033507, 3052618, 542882, 170000, 1030604, 542883, 3052619, 3628000))



# ------------------------------------------------------------------------------------- #

taxID = 550542
name <- withExtend %>% filter(taxid == taxID) %>% pull(binomial) %>% unique()

testHGT <- withExtend %>%
	filter(taxid == taxID & Type == "HGT") %>%
	pull(relGeneStart) %>%
	density(adjust = 1/8, n = 2^12)

testAll <- withExtend %>%
	filter(taxid == taxID) %>%
	pull(relGeneStart) %>%
	density(adjust = 1/8, n = 2^12)


enrichUpColor	<- wes_palette("Darjeeling1")[2]
enrichDownColor	<- wes_palette("Darjeeling1")[1]

plot(testHGT)
plot(testAll)
polys_plus <- lapply(polyclip(A = list("x" = testHGT$x, "y" = testHGT$y), B = list("x" = testAll$x, "y" = testAll$y), op = "minus"), function(poly) {
	tibble(x = poly$x, y = poly$y)
}) %>% bind_rows(.id = "id")

polys_minus <- lapply(polyclip(A = list("x" = testAll$x, "y" = testAll$y), B = list("x" = testHGT$x, "y" = testHGT$y), op = "minus"), function(poly) {
	tibble(x = poly$x, y = poly$y)
}) %>% bind_rows(.id = "id")



MC52_plot <- ggplot() + 
	geom_polygon(data = polys_plus, aes(x = x, y = y, group = id, fill = "Up"), alpha = 0.8) +
	geom_polygon(data = polys_minus, aes(x = x, y = y, group = id, fill = "Down"), alpha = 0.8) +
	geom_vline(aes(xintercept = mauveBounds_mc52$coord/mc_52_max)) +
	geom_label(data = mauveBounds_mc52 %>% filter(orient == "start"), aes(x = coord/mc_52_max, y = 0.2, label = block)) +
	scale_fill_manual(values = c(enrichDownColor, enrichUpColor)) +
	coord_cartesian(xlim = c(0, 1), expand = FALSE) +
	labs(title = name) +
	lightTheme

# ------------------------------------------------------------------------------------- #


taxID = 544556
name <- withExtend %>% filter(taxid == taxID) %>% pull(binomial) %>% unique()

testHGT <- withExtend %>%
	filter(taxid == taxID & Type == "HGT") %>%
	pull(relGeneStart) %>%
	density(adjust = 1/8, n = 2^12)

testAll <- withExtend %>%
	filter(taxid == taxID) %>%
	pull(relGeneStart) %>%
	density(adjust = 1/8, n = 2^12)


enrichUpColor	<- wes_palette("Darjeeling1")[2]
enrichDownColor	<- wes_palette("Darjeeling1")[1]

plot(testHGT)
plot(testAll)
polys_plus <- lapply(polyclip(A = list("x" = testHGT$x, "y" = testHGT$y), B = list("x" = testAll$x, "y" = testAll$y), op = "minus"), function(poly) {
	tibble(x = poly$x, y = poly$y)
}) %>% bind_rows(.id = "id")

polys_minus <- lapply(polyclip(A = list("x" = testAll$x, "y" = testAll$y), B = list("x" = testHGT$x, "y" = testHGT$y), op = "minus"), function(poly) {
	tibble(x = poly$x, y = poly$y)
}) %>% bind_rows(.id = "id")



MC61_plot <- ggplot() + 
	geom_polygon(data = polys_plus, aes(x = x, y = y, group = id, fill = "Up"), alpha = 0.8) +
	geom_polygon(data = polys_minus, aes(x = x, y = y, group = id, fill = "Down"), alpha = 0.8) +
	geom_vline(aes(xintercept = mauveBounds_mc61$coord/mc_61_max)) +
	geom_label(data = mauveBounds_mc61 %>% filter(orient == "start"), aes(x = coord/mc_61_max, y = 0.2, label = block)) +
	scale_fill_manual(values = c(enrichDownColor, enrichUpColor)) +
	coord_cartesian(xlim = c(0, 1), expand = FALSE) +
	labs(title = name) +
	lightTheme

# ------------------------------------------------------------------------------------- #

quartz(width = 12, height = 8, canvas = "white", bg = "white")
print(MC52_plot)
fileName	<- file.path(out_path, "MC52_comp_fig")
invisible(quartz.save(file = paste0(fileName, ".pdf"), type = "pdf", dpi = 100))
invisible(dev.off())

quartz(width = 12, height = 8, canvas = "white", bg = "white")
print(MC61_plot)
fileName	<- file.path(out_path, "MC61_comp_fig")
invisible(quartz.save(file = paste0(fileName, ".pdf"), type = "pdf", dpi = 100))
invisible(dev.off())










