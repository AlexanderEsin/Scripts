invisible(sapply(HGTPos.geo, source, .GlobalEnv))

require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load("tidyverse", "reshape2", "grid", "ggplot2", "ggrepel", "ggdendro", "ggpubr", "wesanderson", "RColorBrewer", "fs")

# ------------------------------------------------------------------------------------- #
message("\nReading in data...", appendLF = FALSE)

# General position data
GPA_perTypeData			<- readRDS(file.path(positionData_path, "AG_perTypeData.rds"))

# DB
dbConn			<- dbConnect(RSQLite::SQLite(), allProtDB_path)

# Zones
GPA_zoneBoundaryList	<- readRDS(file.path(positionData_path, "AG_zoneBoundaries.rds"))

# ------------------------------------------------------------------------------------- #
accass_trans <- read_tsv(file = path(master_path, "Genomes", "Genome_lists", "AG_acc_ass_chromosome.txt"))

phaster_dir <- path(master_path, "PHASTER")
phaster_resDirs <- dir_ls(phaster_dir, type = "directory")


out_path <- "/Users/aesin/Documents/Scripts/HGT_position_reviewer/Figures/phageAnalysis"
dir_create(out_path)

# ------------------------------------------------------------------------------------- #

## Get all the phage detection data together
byGenomePhage_out <- lapply(phaster_resDirs, function(resultDir) {
	
	## Open detail file to get the ID of the genome. First line only
	first_line <- read_lines(file = path(resultDir, "detail.txt"), n_max = 1)
	accession <- unlist(str_split(first_line, "\\|"))[4]
	accAss <- accass_trans %>% filter(Chromosome == accession) %>% pull(AccAss)
	
	# Get the taxid from DB based on genome_l
	getGPA_details	<- dbSendQuery(dbConn, paste0('SELECT DISTINCT taxid, binomial FROM t1 WHERE is_ag == 1 AND acc_ass == "', accAss, '"'))
	GPA_outTbl	<- dbFetch(getGPA_details)
	dbClearResult(getGPA_details)

	message(accAss)
	message(GPA_outTbl$binomial)
	message(resultDir)
	print("\n")
	
	# Get the phage regions by processing the fasta
	fasta_read <- read_lines(file = path(resultDir, "phage_regions.fna"))
	summ_read <- read_lines(file = path(resultDir, "summary.txt"))
	
	pos_lines <- str_extract(str_subset(fasta_read, "^>.+"), "(?<=\\s{2}).+")
	
	# Over each detected phage region, get the start and end positions
	phage_final <- lapply(pos_lines, function(phage_region) {
		split_pos <- unlist(str_split(phage_region, "-"))
		phage_start <- split_pos[[1]]
		phage_end <- split_pos[[2]]
		
		# Get the completeness
		completeness_full <- unlist(str_split(str_subset(summ_read, phage_region), "\\s+"))[4]
		completeness_class <- str_extract(completeness_full, "[a-z]+")
		completeness_score <- str_extract(completeness_full, "[0-9]+")
		
		out_tbl <- GPA_outTbl %>%
			mutate(phage_start = phage_start,
				   phage_end = phage_end,
				   phage_quality = completeness_class,
				   phage_score = completeness_score)
		
		return(out_tbl)
	}) %>% bind_rows() %>%
		mutate(resDir = resultDir,
			   accAss = accAss)
	
	# Get completeness of phage region from summary
	
	
	return(phage_final)
}) %>% bind_rows() %>% as_tibble() %>% type_convert()


# ------------------------------------------------------------------------------------- #

## Get all the protIDs that fall into phage regions per genome - NB some phage regions overlap

GPA_taxids <- byGenomePhage_out %>% distinct(taxid) %>% pull(taxid)

allGPA_phage_data <- lapply(GPA_taxids, function(taxid_id) {
	
	phage_regions <- byGenomePhage_out %>% filter(taxid == taxid_id)
	
	print(taxid_id)
	
	allUniqueInPhage_tbl <- lapply(1:nrow(phage_regions), function(phage_index) {
		
		span <- as_tibble(phage_regions[phage_index,])
		
		gene_start_inRegion <- GPA_perTypeData$All$allPosData %>%
			filter(taxid == span$taxid & geneStart >= span$phage_start & geneStart <= span$phage_end) %>%
			pull(protID)
		
		gene_end_inRegion <- GPA_perTypeData$All$allPosData %>%
			filter(taxid == span$taxid & geneEnd >= span$phage_start & geneEnd <= span$phage_end) %>%
			pull(protID)
		
		genes_partInPhage <- tibble(taxid = taxid_id,
									protID = unique(c(gene_start_inRegion, gene_end_inRegion)),
									phage_start = span$phage_start,
									phage_end = span$phage_end,
									phage_quality = span$phage_quality,
									phage_score = span$phage_score)
		
		return(genes_partInPhage)
	}) %>% bind_rows() %>%
		left_join(phage_regions %>% distinct(taxid, binomial), by = "taxid")
	
	print(head(allUniqueInPhage_tbl))
	

	allUniqueInPhageByType_tbl <- allUniqueInPhage_tbl %>%
		mutate(Type = case_when(
			protID %in% GPA_perTypeData$lHGT$`4`$allPosData$protID ~ "HGT",
			protID %in% GPA_perTypeData$Ver$`3`$allPosData$protID ~ "Ver",
			TRUE ~ "Other"
		))
	
	return(allUniqueInPhageByType_tbl)
}) %>% bind_rows() %>% as_tibble()


# ------------------------------------------------------------------------------------- #

## Genomes with complete phages
GPA_with_confidentPhages <- allGPA_phage_data %>% filter(phage_quality == "intact")

intactPhageGeneCount_plot <- GPA_with_confidentPhages %>%
	mutate(Type = factor(Type, levels = c("HGT", "Ver", "Other"))) %>%
	ggplot() +
	geom_bar(aes(x = Type, y = ..count.., fill = Type), position = position_dodge(), width = 0.5, alpha = 0.8) +
	scale_fill_manual(values = c(dataTypeCols$HGT, dataTypeCols$Ver, dataTypeCols$Other)) +
	lightTheme +
	theme(panel.grid.major.x = element_blank())


## For plotting density with extended coordinates
withRelStart <- GPA_with_confidentPhages %>%
	left_join(GPA_perTypeData$All$allPosData %>% select(protID, geneStart, relGeneStart), by = "protID")

withExtend <- withRelStart %>%
	filter(relGeneStart < 0.2 | relGeneStart > 0.8) %>%
	mutate(relGeneStart = case_when(
		relGeneStart < 0.2 ~ relGeneStart + 1,
		TRUE ~ relGeneStart - 1)
	) %>%
	bind_rows(withRelStart)
	
intactPhageGeneAggregate_GPA_plot <- withExtend %>%
	mutate(Type = factor(Type, levels = c("HGT", "Ver", "Other"))) %>%
	ggplot() +
	geom_rect(data = GPA_zoneBoundaryList$fullRange,
			  aes(xmin = zoneMin, xmax = zoneMax, ymin = -Inf, ymax = Inf),
			  fill = GPA_zoneBoundaryList$fullRange$zoneCol_alpha,
			  inherit.aes = FALSE) +
	stat_density(aes(x = relGeneStart, y = ..count.., fill = Type), colour = "black", geom = "area", size = 0.2, adjust = 1/10, alpha = 0.3, position = "identity") +
	scale_colour_manual(values = c(dataTypeCols$HGT, dataTypeCols$Ver, dataTypeCols$Other)) +
	scale_fill_manual(values = c(dataTypeCols$HGT, dataTypeCols$Ver, dataTypeCols$Other)) +
	coord_cartesian(xlim = c(0, 1), expand = FALSE) +
	lightTheme
	


# ------------------------------------------------------------------------------------- #
## ALL phages

allPhageGeneCount_plot <- allGPA_phage_data %>%
	mutate(Type = factor(Type, levels = c("HGT", "Ver", "Other"))) %>%
	ggplot() +
	geom_bar(aes(x = Type, y = ..count.., fill = Type), position = position_dodge(), width = 0.5, alpha = 0.8) +
	scale_fill_manual(values = c(dataTypeCols$HGT, dataTypeCols$Ver, dataTypeCols$Other)) +
	lightTheme +
	theme(panel.grid.major.x = element_blank())



## For plotting density with extended coordinates
withRelStart <- allGPA_phage_data %>%
	left_join(GPA_perTypeData$All$allPosData %>% select(protID, geneStart, relGeneStart), by = "protID")

withExtend <- withRelStart %>%
	filter(relGeneStart < 0.2 | relGeneStart > 0.8) %>%
	mutate(relGeneStart = case_when(
		relGeneStart < 0.2 ~ relGeneStart + 1,
		TRUE ~ relGeneStart - 1)
	) %>%
	bind_rows(withRelStart)

allPhageGeneAggregate_GPA_plot <- withExtend %>%
	mutate(Type = factor(Type, levels = c("HGT", "Ver", "Other"))) %>%
	ggplot() +
	geom_rect(data = GPA_zoneBoundaryList$fullRange,
			  aes(xmin = zoneMin, xmax = zoneMax, ymin = -Inf, ymax = Inf),
			  fill = GPA_zoneBoundaryList$fullRange$zoneCol_alpha,
			  inherit.aes = FALSE) +
	stat_density(aes(x = relGeneStart, y = ..count.., fill = Type), colour = "black", geom = "area", size = 0.2, adjust = 1/10, alpha = 0.3, position = "identity") +
	scale_colour_manual(values = c(dataTypeCols$HGT, dataTypeCols$Ver, dataTypeCols$Other)) +
	scale_fill_manual(values = c(dataTypeCols$HGT, dataTypeCols$Ver, dataTypeCols$Other)) +
	coord_cartesian(xlim = c(0, 1), expand = FALSE) +
	lightTheme


# ------------------------------------------------------------------------------------- #
# Print outputs

quartz(width = 6, height = 6, canvas = "white", bg = "white")
print(allPhageGeneCount_plot)
fileName	<- file.path(out_path, "intactPhageGeneCount")
invisible(quartz.save(file = paste0(fileName, ".pdf"), type = "pdf", dpi = 100))
invisible(quartz.save(file = paste0(fileName, ".png"), type = "png", dpi = 100))
invisible(dev.off())

quartz(width = 6, height = 6, canvas = "white", bg = "white")
print(allPhageGeneCount_plot)
fileName	<- file.path(out_path, "allPhageGeneCount")
invisible(quartz.save(file = paste0(fileName, ".pdf"), type = "pdf", dpi = 100))
invisible(quartz.save(file = paste0(fileName, ".png"), type = "png", dpi = 100))
invisible(dev.off())


quartz(width = 14, height = 10, canvas = "white", bg = "white")
print(intactPhageGeneAggregate_GPA_plot)
fileName	<- file.path(out_path, "intactPhageAggDistrib")
invisible(quartz.save(file = paste0(fileName, ".pdf"), type = "pdf", dpi = 100))
invisible(quartz.save(file = paste0(fileName, ".png"), type = "png", dpi = 100))
invisible(dev.off())


quartz(width = 14, height = 10, canvas = "white", bg = "white")
print(allPhageGeneAggregate_GPA_plot)
fileName	<- file.path(out_path, "allPhageAggDistrib")
invisible(quartz.save(file = paste0(fileName, ".pdf"), type = "pdf", dpi = 100))
invisible(quartz.save(file = paste0(fileName, ".png"), type = "png", dpi = 100))
invisible(dev.off())


dbDisconnect(dbConn)









