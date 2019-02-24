invisible(sapply(HGTPos.geo, source, .GlobalEnv))

require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load("dplyr", "reshape2", "grid", "ggplot2", "ggrepel", "ggdendro", "ggpubr", "wesanderson", "RColorBrewer", "fs")

bac_rds_path <- "/Users/aesin/Desktop/Bacillus/HGT_position/DataObjects/"
sta_rds_path <- "/Users/aesin/Desktop/Staph/HGT_position/DataObjects/"

out_path <- "/Users/aesin/Documents/Scripts/HGT_position_reviewer/Figures/vertical_gene_synteny"
dir_create(out_path)

# ------------------------------------------------------------------------------------- #
message("\nReading in data...", appendLF = FALSE)

# General position data
GPA_perTypeData			<- readRDS(file.path(positionData_path, "AG_perTypeData.rds"))
BAC_perTypeData			<- readRDS(path(bac_rds_path, "AG_perTypeData.rds"))
STA_perTypeData			<- readRDS(path(sta_rds_path, "AG_perTypeData.rds"))

# Subgroup data
GPA_dnaA_pos_data		<- readRDS(file.path(positionData_path, "bySpecies_dnaA_data.rds"))

# Zone boundaries
GPA_zoneBoundaryList	<- readRDS(file.path(positionData_path, "AG_zoneBoundaries.rds"))

message("\rReading in data... done\n")

# ------------------------------------------------------------------------------------- #

gKau_ver_data <- GPA_perTypeData$Ver$'3'$allPosData %>%
	filter(binomial == "Geobacillus kaustophilus HTA426") %>%
	as_tibble()

bSub_all_data <- BAC_perTypeData$All$allPosData %>%
	filter(binomial == "Bacillus subtilis subsp. subtilis str. 168") %>%
	as_tibble()

sAur_all_data <- STA_perTypeData$All$allPosData %>%
	filter(binomial == "Staphylococcus aureus subsp. aureus NCTC 8325") %>%
	as_tibble()

# ------------------------------------------------------------------------------------- #

dbConn			<- dbConnect(RSQLite::SQLite(), allProtDB_path)

gKau_genomeLen_tbl	<- dbSendQuery(dbConn, 'SELECT DISTINCT genome_l FROM t1 WHERE taxid == 235909')
gKau_genomeLen	<- dbFetch(gKau_genomeLen_tbl) %>% pull(genome_l)
dbClearResult(gKau_genomeLen_tbl)

bSub_allVer_protID_tbl	<- dbSendQuery(dbConn, 'SELECT * FROM t1 WHERE taxid == 224308')
bSub_allGenes_df	<- dbFetch(bSub_allVer_protID_tbl) %>%
	select(-c(sequence, NuclSeq)) %>%
	as_tibble()
dbClearResult(bSub_allVer_protID_tbl)

sAur_allVer_protID_tbl	<- dbSendQuery(dbConn, 'SELECT * FROM t1 WHERE taxid == 93061')
sAur_allGenes_df	<- dbFetch(sAur_allVer_protID_tbl) %>%
	select(-c(sequence, NuclSeq)) %>%
	as_tibble()
dbClearResult(sAur_allVer_protID_tbl)

dbDisconnect(dbConn)

# ------------------------------------------------------------------------------------- #

bSub_gKauVer_uniq <- bSub_allGenes_df %>%
	filter(OrthGroup %in% unique(gKau_ver_data$orthGroup)) %>%
	group_by(OrthGroup) %>%
	summarise(n = n()) %>%
	filter(n == 1)

bSub_gKauVer_data <- bSub_allGenes_df %>%
	filter(OrthGroup %in% bSub_gKauVer_uniq$OrthGroup) %>%
	select(protID, oldOrthGroup = OrthGroup) %>%
	left_join(bSub_all_data, by = "protID") %>%
	select(-orthGroup) %>%
	rename(orthGroup = oldOrthGroup)


sAur_gKauVer_uniq <- sAur_allGenes_df %>%
	filter(OrthGroup %in% unique(gKau_ver_data$orthGroup)) %>%
	group_by(OrthGroup) %>%
	summarise(n = n()) %>%
	filter(n == 1)

sAur_gKauVer_data <- sAur_allGenes_df %>%
	filter(OrthGroup %in% sAur_gKauVer_uniq$OrthGroup) %>%
	select(protID, oldOrthGroup = OrthGroup) %>%
	left_join(sAur_all_data, by = "protID") %>%
	select(-orthGroup) %>%
	rename(orthGroup = oldOrthGroup)

# ------------------------------------------------------------------------------------- #

gKau_in_bSub <- gKau_ver_data %>%
	filter(orthGroup %in% bSub_gKauVer_data$orthGroup)

gKau_in_sAur <- gKau_ver_data %>%
	filter(orthGroup %in% sAur_gKauVer_data$orthGroup)



gKau_bSub_allVer_comb <- bind_rows(gKau_in_bSub, bSub_gKauVer_data) %>% 
	select(-CircStart, -CircEnd) %>%
	mutate(genomeLength = case_when(
		taxid == 235909 ~ gKau_genomeLen,
		TRUE ~ unique(bSub_allGenes_df$genome_l))) %>%
	mutate(relTerStart = geneStart - (genomeLength / 2))

quartz(width = 21, height = 8, canvas = "white", bg = "white")
allVer_bSub_gKau_comparison_plot	<- gKau_bSub_allVer_comb %>%
	ggplot(aes(x = binomial, y = relTerStart, group = orthGroup, color = binomial)) +
	geom_point(shape = 124, size = 8) +
	geom_line(linetype = "solid", color = alpha(wes_palette("IsleofDogs1")[6], 0.9), size = 0.1) +
	scale_color_manual(values = wes_palette("Darjeeling1")[c(5,1)], guide = FALSE) +
	coord_flip() +
	labs(title = "Synteny of vertical genes between Geobacillus kaustophilus and Bacillus subtilis str 168", 
		 subtitle = "Based on Gkau vertical genes conserved in B. subtilis",
		 caption = paste("Number of genes =", nrow(gKau_bSub_allVer_comb) / 2)) +
	theme_classic() +
	theme(
		panel.grid.major.y = element_line(size = 0.05),
		panel.grid.major.x  = element_line(size = 0.4, linetype = "dashed", color = wes_palette("Darjeeling1")[2]),
		axis.text = element_text(size = 11),
		axis.text.y = element_text(angle = 45),
		axis.title = element_text(size = 14),
		axis.title.y = element_blank(),
		plot.title = element_text(hjust = 0.5, size = 16))
print(allVer_bSub_gKau_comparison_plot)

fileName	<- file.path(out_path, "allVer_bSub_gKau_terCentred")
invisible(quartz.save(file = paste0(fileName, ".png"), type = "png", dpi = 100))
invisible(dev.off())



sAur_bSub_allVer_comb <- bind_rows(gKau_in_sAur, sAur_gKauVer_data) %>% 
	select(-CircStart, -CircEnd) %>%
	mutate(genomeLength = case_when(
		taxid == 235909 ~ gKau_genomeLen,
		TRUE ~ unique(sAur_allGenes_df$genome_l))) %>%
	mutate(relTerStart = geneStart - (genomeLength / 2))

quartz(width = 21, height = 8, canvas = "white", bg = "white")
allVer_sAur_gKau_comparison_plot <- sAur_bSub_allVer_comb %>%
	mutate(binomial = factor(binomial, levels = rev(unique(binomial)))) %>%
	ggplot(aes(x = binomial, y = relTerStart, group = orthGroup, color = binomial)) +
	geom_point(shape = 124, size = 8) +
	geom_line(linetype = "solid", color = alpha(wes_palette("IsleofDogs1")[6], 0.9), size = 0.1) +
	scale_color_manual(values = wes_palette("Darjeeling1")[c(5,1)], guide = FALSE) +
	coord_flip() +
	labs(title = "Synteny of vertical genes between Geobacillus kaustophilus and Staph. aureus", 
		 subtitle = "Based on Gkau vertical genes conserved in S. aureus",
		 caption = paste("Number of genes =", nrow(sAur_bSub_allVer_comb) / 2)) +
	theme_classic() +
	theme(
		panel.grid.major.y = element_line(size = 0.05),
		panel.grid.major.x  = element_line(size = 0.4, linetype = "dashed", color = wes_palette("Darjeeling1")[2]),
		axis.text = element_text(size = 11),
		axis.text.y = element_text(angle = 45),
		axis.title = element_text(size = 14),
		axis.title.y = element_blank(),
		plot.title = element_text(hjust = 0.5, size = 16))
print(allVer_sAur_gKau_comparison_plot)

fileName	<- file.path(out_path, "allVer_sAur_gKau_terCentred")
invisible(quartz.save(file = paste0(fileName, ".pdf"), type = "pdf", dpi = 100))
invisible(quartz.save(file = paste0(fileName, ".png"), type = "png", dpi = 100))
invisible(dev.off())




