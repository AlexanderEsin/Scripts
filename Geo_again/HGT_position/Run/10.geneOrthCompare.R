#!/usr/bin/env Rscript

# Load master variables and HGT position functions
invisible(sapply(HGTPos.all, source, .GlobalEnv))

# Any libraries explicitly used in the script
require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load("RSQLite", "dplyr","ggplot2", "ggpubr", "wesanderson", "gridExtra", "parallel", "data.table")


# -------------------------------------- #
# Set global quartz options
quartz.options(canvas = "white", bg = "white")

# Open All_prot database
dbConn			<- dbConnect(RSQLite::SQLite(), allProtDB_path)

# -------------------------------------- #
# Read in data

message("\nReading in data...", appendLF = FALSE)

# Zone list
zoneBoundaryList	<- readRDS(file.path(positionData_path, "AG_zoneBoundaries.rds"))

message("\rReading in data... done\n")

# -------------------------------------- #
# Output path for figures
orthPositionFig_path	<- file.path(figureOutput_path, "orthPosition")
if(!dir.exists(orthPositionFig_path)) dir.create(orthPositionFig_path)

# ------------------------------------------------------------------------------------- #

if (file.exists(file.path(positionData_path, "soloHGTEvent_data.rds"))) {
	soloEventData_df	<- readRDS(file.path(positionData_path, "soloHGTEvent_data.rds"))
} else {

	# List of all the unique COG categorties found in the entire AG dataset
	uniqueCOGs		<- unique(unlist(perTypeData$All$allPosData$COGcat))

	# For each COG and each individual HGT event belonging to that COG, find the range of distance to origin
	# Calculate the max range and SD of the gene positions per-event/per-COG
	byCOG_geneSpread_list	<- mclapply(uniqueCOGs, function(COG) {

		perCOG_data		<- perTypeCOG_data$lHGT$perCOGdata[[COG]]$allData
		eventUnique		<- unique(perCOG_data$eventIndex)
		eventNumber		<- length(eventUnique)

		byEventSpread_list	<- lapply(eventUnique, function(hgtEvent) {

			# Extract just data corresponding to this event
			perEvent_data	<- perCOG_data[which(perCOG_data$eventIndex == hgtEvent),]

			# Number of genes involved in HGT event
			numTotal	<- nrow(perEvent_data)

			# Calculate relative distance to origin (ori-ter symmetry agnostic) and calculate SD on that
			oriDists		<- unlist(lapply(perEvent_data$relGeneStart, function(geneStart) ifelse(geneStart > 0.5, 1 - geneStart, geneStart)))
			eventDistSD		<- sd(perEvent_data$distToOri)
			eventDistMax	<- max(perEvent_data$distToOri) - min(perEvent_data$distToOri)

			# All members of event must belong to single Orth group. Sanity check
			orthGroup	<- unique(perEvent_data$orthGroup)
			if (length(orthGroup) > 1) {
				stop("Should not have more than one group for HGT event")
			}

			# Summary df
			summary_df	<- data.frame(COG = COG, orthGroup = orthGroup, eventIndex = hgtEvent, totalGenes = numTotal, geneSD = eventDistSD, geneDistMax = eventDistMax, stringsAsFactors = FALSE)
			return(list(allData = perEvent_data, summary_df = summary_df))
		})

		# Combine the per-event results (into per-COG results)
		byEventAllData_df 	<- bind_rows(lapply(byEventSpread_list, function(element) return(element$allData)))
		byEventSummary_df 	<- bind_rows(lapply(byEventSpread_list, function(element) return(element$summary_df)))
		return(list(allData = byEventAllData_df, summaryData = byEventSummary_df))
	}, mc.cores = 14)
	# Combine the cross-COG results
	byCOG_geneSpreadAll_df		<- bind_rows(lapply(byCOG_geneSpread_list, function(element) return(element$allData)))
	byCOG_geneSpreadSummary_df	<- bind_rows(lapply(byCOG_geneSpread_list, function(element) return(element$summaryData)))


	# Which COGs have the most independent events?
	eventNumberByCOG	<- sort(table(byCOG_geneSpreadSummary_df$COG), decreasing = TRUE)
	eventNumberByCOG_df	<- data.frame(eventNumberByCOG)
	names(eventNumberByCOG_df)	<- c("COG", "Freq")

	# Sort the COG categories by the number of events
	byCOG_geneSpreadSummary_df$COG	<- factor(byCOG_geneSpreadSummary_df$COG, levels = levels(eventNumberByCOG_df$COG))

	# Check which orthologous groups occur more than once - indicating multiple events in same gene family
	multiEventGroups_list	<- unique(byCOG_geneSpreadSummary_df$orthGroup[duplicated(byCOG_geneSpreadSummary_df$orthGroup)])

	# Entries where there is only a single HGT event per orthGroup
	soloEventGroups_list	<- byCOG_geneSpreadSummary_df[which(!byCOG_geneSpreadSummary_df$orthGroup %in% multiEventGroups_list),]
	soloEventData_df		<- byCOG_geneSpreadAll_df[which(byCOG_geneSpreadAll_df$orthGroup %in% soloEventGroups_list$orthGroup),]

	saveRDS(soloEventData_df, file = file.path(positionData_path, "soloHGTEvent_data.rds"))

}

# ------------------------------------------------------------------------------------- #


dnaA_clean_file	<- file.path(positionData_path, "bySpecies_dnaA_data.rds")
if (file.exists(dnaA_clean_file)) {
	dnaA_clean_trim	<- readRDS(dnaA_clean_file)

	# Read in the species subgroupings file (result of Scripts/Geo_again/Genomes/8.Identify_species_subgroupings.R)
	speciesSubgroups_file	<- file.path(genome_path, "Genome_lists", "Species_groupings.tsv")
	speciesSubgroups_df		<- read.table(speciesSubgroups_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "")

} else {
	# Extract all the dnaA genes for the 
	dnaA_protID_list_file	<- file.path(position_path, "dnaA_protIDs.txt")
	dnaA_protIDs			<- read.table(file = dnaA_protID_list_file, sep = "\n", header = FALSE, stringsAsFactors = FALSE)$V1

	# 4982 entries in this table, of which 4975 are unique (so 7 genomes have 2x dnaA gene?)
	dnaA_data		<- dbSendQuery(dbConn, 'SELECT * FROM t1 WHERE protID = :protIDs')
	dbBind(dnaA_data, param = list(protIDs = dnaA_protIDs))
	dnaA_data_df	<- dbFetch(dnaA_data)
	dbClearResult(dnaA_data)

	# Clean DNAa table to contain only taxid entries with one DNAa gene
	duplDNAa_taxids	<- dnaA_data_df$taxid[which(duplicated(dnaA_data_df$taxid))]
	dnaA_clean_df	<- dnaA_data_df[which(!dnaA_data_df$taxid %in% duplDNAa_taxids),]

	# Clean dnaA table to not contain entries that are designated to be on plasmids
	dnaA_clean_df	<- dnaA_clean_df[which(!dnaA_clean_df$plasmid == "T"),]

	# Read in the species subgroupings file (result of Scripts/Geo_again/Genomes/8.Identify_species_subgroupings.R)
	speciesSubgroups_file	<- file.path(genome_path, "Genome_lists", "Species_groupings.tsv")
	speciesSubgroups_df		<- read.table(speciesSubgroups_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "")

	# Select just the columns we need
	dnaA_clean_trim	<- subset(dnaA_clean_df, select = c(taxid, gene_start, gene_end, strand))
	names(dnaA_clean_trim)	<- c("taxid", "oriStart", "oriEnd", "oriStrand")

	saveRDS(object = dnaA_clean_trim, file = file.path(positionData_path, "bySpecies_dnaA_data.rds"))

	gc()
}

# ------------------------------------------------------------------------------------- #
# "Bacillaceae"
sporeGroups		<- c("Bacillaceae", "Paenibacillaceae", "Alicyclobacillaceae", "Pasteuriaceae", "Sporolactobacillaceae", "Thermoactinomycetaceae")
nonSporeGroups	<- c("Listeriaceae", "Staphylococcaceae", "Lactobacillales")

# Check which groups have not been done
allGroups		<- c(sporeGroups, nonSporeGroups)
groupsToPrepare	<- lapply(allGroups, function(taxonGroup) {
	file_name	<- paste0(taxonGroup, ".rds")
	file_path	<- file.path(supergroupData_path, file_name)
	if (!file.exists(file_path)) {
		return(taxonGroup)
	} else {
		return(NA)
	}
})
groupsToPrepare	<- groupsToPrepare[!is.na(groupsToPrepare)]

# If we need to prepare groups, read in the Nodes and Names data tables
if (length(groupsToPrepare) > 0) {
	nodes_dt	<- readRDS(file.path(taxdmp_path, "nodes.rds"))
	names_dt	<- readRDS(file.path(taxdmp_path, "names.rds"))

	# Prepare and save the groups that haven't been done before
	invisible(lapply(groupsToPrepare, function(taxonGroup) {

		# Identify the parentSpecies taxids that belong to desired superGroups
		superGroup_taxids	<- taxidsBySupergroup(taxonGroup, speciesGroupings = speciesSubgroups_df, nodes_dt = nodes_dt, names_dt = names_dt)

		# Save the result
		outputName	<- file_name	<- paste0(taxonGroup, ".rds")
		saveRDS(object = superGroup_taxids, file = file.path(supergroupData_path, outputName))

	}))

	# Clean away the names and nodes data tables
	rm(list = c("names_dt", "nodes_dt")); gc()

}


superGroup_readList	<- lapply(allGroups, function(taxonGroup) {
	file_name	<- paste0(taxonGroup, ".rds")
	file_path	<- file.path(supergroupData_path, file_name)

	taxonData	<- readRDS(file_path)
	return(taxonData)
})
names(superGroup_readList)	<- allGroups


# ------------------------------------------------------------------------------------- #

allSoloCOGCat	<- unlist(unique(soloEventData_df$COGcat))

# All ortholog data
perCOGOrth_data			<- lapply(allSoloCOGCat, processOrthologPosition,
	orthData = soloEventData_df,
	dbConn = dbConn,
	dnaA_data = dnaA_clean_trim,
	speciesSubgroups = speciesSubgroups_df,
	superGroupTaxid_list = superGroup_readList,
	outlierTaxid = outlierTaxid)
names(perCOGOrth_data)	<- allSoloCOGCat
gc()


# ------------------------------------------------------------------------------------- #

all_posBias_dt	<- as.data.table(bind_rows(lapply(perCOGOrth_data, function(COG) return(COG$supergroupCompareData))))
all_posBias_dt	<- all_posBias_dt %>% 
	group_by(COGcat) %>%
	mutate(sporeForming = case_when(
			SuperGroup %in% sporeGroups ~ "Spore",
			SuperGroup %in% nonSporeGroups ~ "No Spore",
			SuperGroup == "HGT" ~ "HGT"
			)
	) %>% 
	mutate(sporeForming = factor(sporeForming, levels = c("Spore", "No Spore", "HGT")))

 # Spore No Spore      HGT 
 #   35172    18883     6946 

# -------------------------------------- #

sporeCols	<- c(wes_palette("Chevalier1")[1:2], wes_palette("FantasticFox1")[5])

cutOffCogs				<- all_posBias_dt %>% group_by(COGcat) %>% filter(sporeForming == "HGT") %>% summarise(count = n()) %>% filter(count >= 100) %>% pull(COGcat)
all_posBias_numCut_dt	<- all_posBias_dt %>% filter(COGcat %in% cutOffCogs)

all_posBias_numCut_ext	<- all_posBias_numCut_dt %>%
	subset(distToOri < 0.1 | distToOri > 0.4) %>%
	mutate(distToOri = case_when(
		distToOri < 0.1 ~ distToOri + 0.5,
		TRUE ~ distToOri - 0.5)
	) %>%
	bind_rows(all_posBias_numCut_dt)


allCOGs_sporeCompare_plot	<- ggplot(data = all_posBias_numCut_ext %>% filter(sporeForming != "HGT"), mapping = aes(x = distToOri, y = ..count.., color = sporeForming, fill = sporeForming)) +
	facet_wrap(~COGcat, scales = "free_y") +
	scale_y_continuous(name = "Gene density") +
	scale_x_continuous(name = "Relative distance from origin") +
	stat_density(geom = "line", adjust = 1/10, size = 0.7, position = "identity") +
	scale_color_manual(values = sporeCols) +
	scale_fill_manual(values = alpha(sporeCols, 0), guide = FALSE) +
	geom_rect(
		data = zoneBoundaryList$halfGenomeRange,
		aes(xmin = zoneMin, xmax = zoneMax, ymin = -Inf, ymax = Inf),
		fill = rep(zoneBoundaryList$halfGenomeRange$zoneCol_alpha, length(cutOffCogs)),
		inherit.aes = FALSE) +
	guides(colour = guide_legend(
		title = "Gene Group",
		label.position = "right",
		title.position = "top")) +
	coord_cartesian(xlim = c(0, 0.5), expand = FALSE) +
	lightTheme + 
	theme(
		panel.spacing = unit(2, "lines"),
		panel.grid.minor = element_blank(),
		legend.position = "right",
		legend.background = element_rect(color = axisCol, fill = "white"),
		legend.text.align = 0.5,
		legend.text = element_text(margin = margin(t = 5, b = 5)),
		strip.background = element_rect(color = axisCol, fill = "white"),
		strip.text = element_text(color = textCol, size = 12)
	)


onlyM	<- all_posBias_numCut_ext %>%
	filter(COGcat == "M")

inSpore_uniqueOrth	<- onlyM %>% filter(sporeForming == "Spore") %>% pull(Group) %>% unique()
nonSpore_uniqueOrth	<- onlyM %>% filter(sporeForming == "No Spore") %>% pull(Group) %>% unique()

onlyInSpore	<- inSpore_uniqueOrth[-which(inSpore_uniqueOrth %in% nonSpore_uniqueOrth)]


quartz(width = 24, height = 12)
print(allCOGs_sporeCompare_plot)
quartz.save(file = file.path(orthPositionFig_path, "allCOG_sporeVSnonspore_byCount.pdf"), type = "pdf", dpi = 300)
invisible(dev.off())

# -------------------------------------- #

COG_M_subset	<- all_posBias_dt %>% filter(COGcat == "M")
COG_M_count		<- COG_M_subset %>% group_by(sporeForming) %>% summarise(nSpecies = n(), nGenes = sum(childNumber)) %>% mutate(label = paste0(sporeForming, "\nGenes = ", nGenes))

COG_M_orthologs_sporeCompare_plot	<- ggplot(data = COG_M_subset, mapping = aes(x = distToOri, color = sporeForming, fill = sporeForming)) +
	scale_y_continuous(name = "Gene density") +
	scale_x_continuous(name = "Relative distance from origin") +
	geom_density(adjust = 1/5, size = 0.7, alpha = 0.3, position = "identity") +
	scale_color_manual(values = c(sporeCols), guide = FALSE) +
	scale_fill_manual(values = c(sporeCols), labels = COG_M_count$label) +
	guides(fill = guide_legend(
		title = "Gene Group",
		label.position = "right",
		title.position = "top")) +
	lightTheme + 
	theme(
		panel.spacing = unit(2, "lines"),
		panel.grid.minor = element_blank(),
		legend.position = "right",
		legend.background = element_rect(color = axisCol, fill = "white"),
		legend.text.align = 0.5,
		legend.text = element_text(margin = margin(t = 5, b = 5)),
		strip.background = element_rect(color = axisCol, fill = "white"),
		strip.text = element_text(color = textCol, size = 12)
	)


quartz(width = 18, height = 8)
print(COG_M_orthologs_sporeCompare_plot)
quartz.save(file = file.path(orthPositionFig_path, "COG_M_sporeVSnonspore.pdf"), type = "pdf", dpi = 300)
invisible(dev.off())




















# -------------------------------------- #
# General position data
perTypeData			<- readRDS(file.path(positionData_path, "AG_perTypeData.rds"))
bacillacInclude		<- FALSE

superGroup_list <- superGroup_readList
if (!bacillacInclude) {
	superGroup_list$Bacillaceae <- NULL
	superGroup_list$Sporolactobacillaceae <- NULL
	superGroup_list$Thermoactinomycetaceae <- NULL
}



# -------------------------------------- #
if (length(superGroup_list[["Bacillaceae"]]) != 0) {
	orthPositionFig_extPath	<- file.path(orthPositionFig_path, "Bacillaceae_included")
} else {
	orthPositionFig_extPath	<- file.path(orthPositionFig_path, "Bacillaceae_excluded")
}

if (!dir.exists(orthPositionFig_extPath)) dir.create(orthPositionFig_extPath)



# -------------------------------------- #



# Isolate all GPA COG M genes
all_GPA_cogM_genes	<- perTypeData$All$allPosData %>% filter(COGcat == "M") %>% mutate(type = "All")
hgt_GPA_cogM_genes	<- perTypeData$lHGT$'4'$allPosData %>% filter(COGcat == "M") %>% mutate(type = "HGT")
ver_GPA_cogM_genes	<- perTypeData$Ver$'3'$allPosData %>% filter(COGcat == "M") %>% mutate(type = "Ver")

GPA_cogM_byType		<- bind_rows(list(all_GPA_cogM_genes, hgt_GPA_cogM_genes, ver_GPA_cogM_genes))


# Expand for boundary effects in plotting
GPA_cogM_byType_ext	<- GPA_cogM_byType %>%
	subset(distToOri < 0.1 | distToOri > 0.4) %>%
	mutate(distToOri = case_when(distToOri < 0.1 ~ distToOri + 0.5, TRUE ~ distToOri - 0.5)) %>%
	bind_rows(GPA_cogM_byType)

# Plot the relationship
all_GPA_cogM_position_plot	<- ggplot(data = GPA_cogM_byType_ext, mapping = aes(x = distToOri, color = type)) +
	scale_y_continuous(name = "Gene density") +
	scale_x_continuous(name = "Relative distance from origin") +
	stat_density(geom = "line", adjust = 1/8, size = 2, alpha = 1, position = "identity") +
	scale_color_manual(values = c(dataTypeCols$All, dataTypeCols$HGT, dataTypeCols$Ver)) +
	geom_rect(
		data = zoneBoundaryList$halfGenomeRange,
		aes(xmin = zoneMin, xmax = zoneMax, ymin = -Inf, ymax = Inf),
		fill = zoneBoundaryList$halfGenomeRange$zoneCol_alpha,
		inherit.aes = FALSE) +
	coord_cartesian(xlim = c(0, 0.5), expand = FALSE) +
	lightTheme

quartz(width = 18, height = 8)
print(all_GPA_cogM_position_plot)
quartz.save(file = file.path(orthPositionFig_path, "COG_M_density_GPA.pdf"), type = "pdf", dpi = 300)
invisible(dev.off())

# Deconstruct plot to find nadir
deconstructPlot <- ggplot_build(all_GPA_cogM_position_plot)
allData_only	<- deconstructPlot$data[[1]] %>% filter(group == 3)
nadirDist		<- allData_only[which(allData_only$scaled == min(allData_only$scaled)),] %>% pull(x)

# -------------------------------------- #
# Select all the M genes in the potential forespore region
firstThird_genes	<- GPA_cogM_byType %>% filter(distToOri < nadirDist & type == "All")

firstThird_orthologs <- processOrthologPosition_2(
	COGname = "M",
	orthData = firstThird_genes,
	dbConn = dbConn,
	dnaA_data = dnaA_clean_trim,
	speciesSubgroups = speciesSubgroups_df,
	superGroupTaxid_list = superGroup_list,
	outlierTaxid = outlierTaxid)


firstThird_posBias_dt	<- as.data.table(firstThird_orthologs$supergroupCompareData)
firstThird_posBias_dt %<>% 
	mutate(sporeForming = case_when(
			SuperGroup %in% sporeGroups ~ "Spore",
			SuperGroup %in% nonSporeGroups ~ "No Spore",
			TRUE ~ "NA")
	) %>% 
	mutate(sporeForming = factor(sporeForming, levels = c("Spore", "No Spore")))


# Select those entries where Orths are found in both Spore and Non-spore and those where only in one
firstThird_posBias_dt %<>%
	group_by(Group) %>% 
	mutate(sporePresence = case_when(
		n_distinct(sporeForming) == 2 ~ "Both",
		sporeForming == "Spore" ~ "S Only",
		TRUE ~ "NS Only"
		)
	)


# There are 32 / 141 groups that have orthologs in ONLY either spore (29) or non-spore (3) - if Bacillaceae are excluded. Otherwise only 38 S-only groups and no NS-only groups
firstThird_posBias_dt %>% group_by(Group, sporePresence) %>% summarise(n = n()) %>% filter(sporePresence != "Both") %>% group_by(sporePresence) %>% summarise(nrow = n())

# Extend the dataset for plotting
firstThird_posBias_dt_ext	<-  firstThird_posBias_dt %>%
	subset(distToOri < 0.1 | distToOri > 0.4) %>%
	mutate(distToOri = case_when(
		distToOri < 0.1 ~ distToOri + 0.5,
		TRUE ~ distToOri - 0.5)
	) %>%
	bind_rows(firstThird_posBias_dt)


# Plot those genes that have ortgologs in either Spore or Non-spore species but not both
speciesRestrcited_firstThird <- firstThird_posBias_dt_ext %>% group_by(Group, sporePresence) %>% filter(sporePresence != "Both")
speciesRestrcited_firstThird_plot <- ggplot(data = speciesRestrcited_firstThird, mapping = aes(x = distToOri, y = ..count.., color = sporePresence, fill = sporePresence)) +
	scale_y_continuous(name = "Gene density") +
	scale_x_continuous(name = "Relative distance from origin") +
	stat_density(geom = "line", adjust = 1/10, size = 2, position = "identity") +
	scale_color_manual(values = sporeCols) +
	scale_fill_manual(values = alpha(sporeCols, 0), guide = FALSE) +
	geom_rect(
		data = zoneBoundaryList$halfGenomeRange,
		aes(xmin = zoneMin, xmax = zoneMax, ymin = -Inf, ymax = Inf),
		fill = zoneBoundaryList$halfGenomeRange$zoneCol_alpha,
		inherit.aes = FALSE) +
	guides(colour = guide_legend(
		title = "Gene Group",
		label.position = "right",
		title.position = "top")) +
	coord_cartesian(xlim = c(0, 0.5), expand = FALSE) +
	lightTheme + 
	theme(
		panel.spacing = unit(2, "lines"),
		panel.grid.minor = element_blank(),
		legend.position = "right",
		legend.background = element_rect(color = axisCol, fill = "white"),
		legend.text.align = 0.5,
		legend.text = element_text(margin = margin(t = 5, b = 5)),
		strip.background = element_rect(color = axisCol, fill = "white"),
		strip.text = element_text(color = textCol, size = 12)
	)

quartz(width = 18, height = 8)
print(speciesRestrcited_firstThird_plot)
quartz.save(file = file.path(orthPositionFig_extPath, "M_FirstThird_speciesRestricted_density.pdf"), type = "pdf", dpi = 300)
invisible(dev.off())



# Plot all orthologs
firstThird_orth_plot	<- ggplot(data = firstThird_posBias_dt_ext, mapping = aes(x = distToOri, y = ..count.., color = sporeForming, fill = sporeForming)) +
	scale_y_continuous(name = "Gene density") +
	# expand_limits(y = 0) +
	scale_x_continuous(name = "Relative distance from origin") +
	stat_density(geom = "line", adjust = 1/10, size = 2, position = "identity") +
	scale_color_manual(values = sporeCols) +
	scale_fill_manual(values = alpha(sporeCols, 0), guide = FALSE) +
	geom_rect(
		data = zoneBoundaryList$halfGenomeRange,
		aes(xmin = zoneMin, xmax = zoneMax, ymin = -Inf, ymax = Inf),
		fill = zoneBoundaryList$halfGenomeRange$zoneCol_alpha,
		inherit.aes = FALSE) +
	guides(colour = guide_legend(
		title = "Gene Group",
		label.position = "right",
		title.position = "top")) +
	coord_cartesian(xlim = c(0, 0.5), expand = FALSE) +
	lightTheme + 
	theme(
		panel.spacing = unit(2, "lines"),
		panel.grid.minor = element_blank(),
		legend.position = "right",
		legend.background = element_rect(color = axisCol, fill = "white"),
		legend.text.align = 0.5,
		legend.text = element_text(margin = margin(t = 5, b = 5)),
		strip.background = element_rect(color = axisCol, fill = "white"),
		strip.text = element_text(color = textCol, size = 12)
	)

quartz(width = 18, height = 8)
print(firstThird_orth_plot)
quartz.save(file = file.path(orthPositionFig_extPath, "M_All_firstThird_density.pdf"), type = "pdf", dpi = 300)
invisible(dev.off())

# -------------------------------------- #
# Other two thirds
twoThird_genes	<- GPA_cogM_byType %>% filter(distToOri > nadirDist & type == "All")


twoThird_orthologs <- processOrthologPosition_2(
	COGname = "M",
	orthData = twoThird_genes,
	dbConn = dbConn,
	dnaA_data = dnaA_clean_trim,
	speciesSubgroups = speciesSubgroups_df,
	superGroupTaxid_list = superGroup_list,
	outlierTaxid = outlierTaxid)


twoThird_posBias_dt	<- as.data.table(twoThird_orthologs$supergroupCompareData)
twoThird_posBias_dt	<- twoThird_posBias_dt %>% 
	mutate(sporeForming = case_when(
			SuperGroup %in% sporeGroups ~ "Spore",
			SuperGroup %in% nonSporeGroups ~ "No Spore",
			TRUE ~ "NA")
	) %>% 
	mutate(sporeForming = factor(sporeForming, levels = c("Spore", "No Spore")))

# Select those entries where Orths are found in both Spore and Non-spore and those where only in one
twoThird_posBias_dt %<>%
	group_by(Group) %>% 
	mutate(sporePresence = case_when(
		n_distinct(sporeForming) == 2 ~ "Both",
		sporeForming == "Spore" ~ "S Only",
		TRUE ~ "NS Only"
		)
	)

# Extend the dataset for plotting
twoThird_posBias_dt_ext	<-  twoThird_posBias_dt %>%
	subset(distToOri < 0.1 | distToOri > 0.4) %>%
	mutate(distToOri = case_when(
		distToOri < 0.1 ~ distToOri + 0.5,
		TRUE ~ distToOri - 0.5)
	) %>%
	bind_rows(twoThird_posBias_dt)


# Plot density of genes in groups that have Orths restricted to either S or non-S only
speciesRestrcited_twoThird <- twoThird_posBias_dt_ext %>% group_by(Group, sporePresence) %>% filter(sporePresence != "Both")
speciesRestrcited_twoThird_plot <- ggplot(data = speciesRestrcited_twoThird, mapping = aes(x = distToOri, y = ..count.., color = sporePresence, fill = sporePresence)) +
	scale_y_continuous(name = "Gene density") +
	scale_x_continuous(name = "Relative distance from origin") +
	stat_density(geom = "line", adjust = 1/10, size = 2, position = "identity") +
	scale_color_manual(values = sporeCols) +
	scale_fill_manual(values = alpha(sporeCols, 0), guide = FALSE) +
	geom_rect(
		data = zoneBoundaryList$halfGenomeRange,
		aes(xmin = zoneMin, xmax = zoneMax, ymin = -Inf, ymax = Inf),
		fill = zoneBoundaryList$halfGenomeRange$zoneCol_alpha,
		inherit.aes = FALSE) +
	guides(colour = guide_legend(
		title = "Gene Group",
		label.position = "right",
		title.position = "top")) +
	coord_cartesian(xlim = c(0, 0.5), expand = FALSE) +
	lightTheme + 
	theme(
		panel.spacing = unit(2, "lines"),
		panel.grid.minor = element_blank(),
		legend.position = "right",
		legend.background = element_rect(color = axisCol, fill = "white"),
		legend.text.align = 0.5,
		legend.text = element_text(margin = margin(t = 5, b = 5)),
		strip.background = element_rect(color = axisCol, fill = "white"),
		strip.text = element_text(color = textCol, size = 12)
	)

quartz(width = 18, height = 8)
print(speciesRestrcited_twoThird_plot)
quartz.save(file = file.path(orthPositionFig_extPath, "M_TwoThird_speciesRestricted_density.pdf"), type = "pdf", dpi = 300)
invisible(dev.off())


twoThird_orth_plot	<- ggplot(data = twoThird_posBias_dt_ext, mapping = aes(x = distToOri, y = ..count.., color = sporeForming, fill = sporeForming)) +
	scale_y_continuous(name = "Gene density") +
	scale_x_continuous(name = "Relative distance from origin") +
	stat_density(geom = "line", adjust = 1/10, size = 2, position = "identity") +
	scale_color_manual(values = sporeCols) +
	scale_fill_manual(values = alpha(sporeCols, 0), guide = FALSE) +
	geom_rect(
		data = zoneBoundaryList$halfGenomeRange,
		aes(xmin = zoneMin, xmax = zoneMax, ymin = -Inf, ymax = Inf),
		fill = zoneBoundaryList$halfGenomeRange$zoneCol_alpha,
		inherit.aes = FALSE) +
	guides(colour = guide_legend(
		title = "Gene Group",
		label.position = "right",
		title.position = "top")) +
	coord_cartesian(xlim = c(0, 0.5), expand = FALSE) +
	lightTheme + 
	theme(
		panel.spacing = unit(2, "lines"),
		panel.grid.minor = element_blank(),
		legend.position = "right",
		legend.background = element_rect(color = axisCol, fill = "white"),
		legend.text.align = 0.5,
		legend.text = element_text(margin = margin(t = 5, b = 5)),
		strip.background = element_rect(color = axisCol, fill = "white"),
		strip.text = element_text(color = textCol, size = 12)
	)
quartz(width = 18, height = 8)
print(twoThird_orth_plot)
quartz.save(file = file.path(orthPositionFig_extPath, "M_All_twoThird_density.pdf"), type = "pdf", dpi = 300)
invisible(dev.off())




# -------------------------------------- #
# All combined

allMOrths_pos	<- bind_rows(firstThird_posBias_dt, twoThird_posBias_dt)

allMOrths_pos %>% group_by(Group, sporePresence) %>% summarise(n = n()) %>% filter(sporePresence != "Both") %>% group_by(sporePresence) %>% summarise(nrow = n())


# Extend the dataset for plotting
allMOrths_pos_ext	<-  allMOrths_pos %>%
	subset(distToOri < 0.1 | distToOri > 0.4) %>%
	mutate(distToOri = case_when(
		distToOri < 0.1 ~ distToOri + 0.5,
		TRUE ~ distToOri - 0.5)
	) %>%
	bind_rows(allMOrths_pos)


allM_orth_plot	<- ggplot(data = allMOrths_pos_ext, mapping = aes(x = distToOri, y = ..count.., color = sporeForming, fill = sporeForming)) +
	scale_y_continuous(name = "Gene density") +
	scale_x_continuous(name = "Relative distance from origin") +
	stat_density(geom = "line", adjust = 1/10, size = 2, position = "identity") +
	scale_color_manual(values = sporeCols) +
	scale_fill_manual(values = alpha(sporeCols, 0), guide = FALSE) +
	geom_rect(
		data = zoneBoundaryList$halfGenomeRange,
		aes(xmin = zoneMin, xmax = zoneMax, ymin = -Inf, ymax = Inf),
		fill = zoneBoundaryList$halfGenomeRange$zoneCol_alpha,
		inherit.aes = FALSE) +
	guides(colour = guide_legend(
		title = "Gene Group",
		label.position = "right",
		title.position = "top")) +
	coord_cartesian(xlim = c(0, 0.5), expand = FALSE) +
	lightTheme + 
	theme(
		panel.spacing = unit(2, "lines"),
		panel.grid.minor = element_blank(),
		legend.position = "right",
		legend.background = element_rect(color = axisCol, fill = "white"),
		legend.text.align = 0.5,
		legend.text = element_text(margin = margin(t = 5, b = 5)),
		strip.background = element_rect(color = axisCol, fill = "white"),
		strip.text = element_text(color = textCol, size = 12)
	)
quartz(width = 18, height = 8)
print(allM_orth_plot)
quartz.save(file = file.path(orthPositionFig_extPath, "Combined_M_density.pdf"), type = "pdf", dpi = 300)
invisible(dev.off())























# -------------------------------------- #
# -------------------------------------- #
# -------------------------------------- #
# -------------------------------------- #























testCOG	= "P"

# Isolate all GPA COG M genes
all_GPA_testCOG_genes	<- perTypeData$All$allPosData %>% filter(COGcat == testCOG) %>% mutate(type = "All")
hgt_GPA_testCOG_genes	<- perTypeData$lHGT$'4'$allPosData %>% filter(COGcat == testCOG) %>% mutate(type = "HGT")
ver_GPA_testCOG_genes	<- perTypeData$Ver$'3'$allPosData %>% filter(COGcat == testCOG) %>% mutate(type = "Ver")

GPA_testCOG_byType		<- bind_rows(list(all_GPA_testCOG_genes, hgt_GPA_testCOG_genes, ver_GPA_testCOG_genes))


# Expand for boundary effects in plotting
GPA_testCOG_byType_ext	<- GPA_testCOG_byType %>%
	subset(distToOri < 0.1 | distToOri > 0.4) %>%
	mutate(distToOri = case_when(distToOri < 0.1 ~ distToOri + 0.5, TRUE ~ distToOri - 0.5)) %>%
	bind_rows(GPA_testCOG_byType)

# Plot the relationship
all_GPA_testCOG_position_plot	<- ggplot(data = GPA_testCOG_byType_ext, mapping = aes(x = distToOri, color = type)) +
	scale_y_continuous(name = "Gene density") +
	scale_x_continuous(name = "Relative distance from origin") +
	stat_density(geom = "line", adjust = 1/8, size = 2, alpha = 1, position = "identity") +
	scale_color_manual(values = c(dataTypeCols$All, dataTypeCols$HGT, dataTypeCols$Ver)) +
	geom_rect(
		data = zoneBoundaryList$halfGenomeRange,
		aes(xmin = zoneMin, xmax = zoneMax, ymin = -Inf, ymax = Inf),
		fill = zoneBoundaryList$halfGenomeRange$zoneCol_alpha,
		inherit.aes = FALSE) +
	coord_cartesian(xlim = c(0, 0.5), expand = FALSE) +
	lightTheme

quartz(width = 18, height = 8)
print(all_GPA_testCOG_position_plot)
quartz.save(file = file.path(orthPositionFig_path, paste0("COG_", testCOG, "_density_GPA.pdf")), type = "pdf", dpi = 300)
invisible(dev.off())


# We use the "Nadir Distance" from the M calculations above

# -------------------------------------- #
# Select all the M genes in the potential forespore region
firstThird_genes	<- GPA_testCOG_byType %>% filter(distToOri < nadirDist & type == "All")

firstThird_orthologs <- processOrthologPosition_2(
	COGname = testCOG,
	orthData = firstThird_genes,
	dbConn = dbConn,
	dnaA_data = dnaA_clean_trim,
	speciesSubgroups = speciesSubgroups_df,
	superGroupTaxid_list = superGroup_list,
	outlierTaxid = outlierTaxid)


firstThird_posBias_dt	<- as.data.table(firstThird_orthologs$supergroupCompareData)
firstThird_posBias_dt	<- firstThird_posBias_dt %>% 
	mutate(sporeForming = case_when(
			SuperGroup %in% sporeGroups ~ "Spore",
			SuperGroup %in% nonSporeGroups ~ "No Spore",
			TRUE ~ "NA")
	) %>% 
	mutate(sporeForming = factor(sporeForming, levels = c("Spore", "No Spore")))


# Select those entries where Orths are found in both Spore and Non-spore and those where only in one
firstThird_posBias_dt %<>%
	group_by(Group) %>% 
	mutate(sporePresence = case_when(
		n_distinct(sporeForming) == 2 ~ "Both",
		sporeForming == "Spore" ~ "S Only",
		TRUE ~ "NS Only"
		)
	)


# There are 32 / 141 groups that have orthologs in ONLY either spore (29) or non-spore (3) - if Bacillaceae are excluded. Otherwise only 38 S-only groups and no NS-only groups
firstThird_posBias_dt %>% group_by(Group, sporePresence) %>% summarise(n = n()) %>% filter(sporePresence != "Both") %>% group_by(sporePresence) %>% summarise(nrow = n())


# Extend the dataset for plotting
firstThird_posBias_dt_ext	<-  firstThird_posBias_dt %>%
	subset(distToOri < 0.1 | distToOri > 0.4) %>%
	mutate(distToOri = case_when(
		distToOri < 0.1 ~ distToOri + 0.5,
		TRUE ~ distToOri - 0.5)
	) %>%
	bind_rows(firstThird_posBias_dt)


# Plot those genes that have ortgologs in either Spore or Non-spore species but not both
speciesRestrcited_firstThird <- firstThird_posBias_dt_ext %>% group_by(Group, sporePresence) %>% filter(sporePresence != "Both")
speciesRestrcited_firstThird_plot <- ggplot(data = speciesRestrcited_firstThird, mapping = aes(x = distToOri, y = ..count.., color = sporePresence, fill = sporePresence)) +
	scale_y_continuous(name = "Gene density") +
	scale_x_continuous(name = "Relative distance from origin") +
	stat_density(geom = "line", adjust = 1/10, size = 2, position = "identity") +
	scale_color_manual(values = sporeCols) +
	scale_fill_manual(values = alpha(sporeCols, 0), guide = FALSE) +
	geom_rect(
		data = zoneBoundaryList$halfGenomeRange,
		aes(xmin = zoneMin, xmax = zoneMax, ymin = -Inf, ymax = Inf),
		fill = zoneBoundaryList$halfGenomeRange$zoneCol_alpha,
		inherit.aes = FALSE) +
	guides(colour = guide_legend(
		title = "Gene Group",
		label.position = "right",
		title.position = "top")) +
	coord_cartesian(xlim = c(0, 0.5), expand = FALSE) +
	lightTheme + 
	theme(
		panel.spacing = unit(2, "lines"),
		panel.grid.minor = element_blank(),
		legend.position = "right",
		legend.background = element_rect(color = axisCol, fill = "white"),
		legend.text.align = 0.5,
		legend.text = element_text(margin = margin(t = 5, b = 5)),
		strip.background = element_rect(color = axisCol, fill = "white"),
		strip.text = element_text(color = textCol, size = 12)
	)

# quartz(width = 18, height = 8)
# print(speciesRestrcited_firstThird_plot)
# quartz.save(file = file.path(orthPositionFig_extPath, paste0(testCOG, "_FirstThird_speciesRestricted_density.pdf")), type = "pdf", dpi = 300)
# invisible(dev.off())



# Plot all orthologs
firstThird_orth_plot	<- ggplot(data = firstThird_posBias_dt_ext, mapping = aes(x = distToOri, y = ..count.., color = sporeForming, fill = sporeForming)) +
	scale_y_continuous(name = "Gene density") +
	scale_x_continuous(name = "Relative distance from origin") +
	stat_density(geom = "line", adjust = 1/10, size = 2, position = "identity") +
	scale_color_manual(values = sporeCols) +
	scale_fill_manual(values = alpha(sporeCols, 0), guide = FALSE) +
	geom_rect(
		data = zoneBoundaryList$halfGenomeRange,
		aes(xmin = zoneMin, xmax = zoneMax, ymin = -Inf, ymax = Inf),
		fill = zoneBoundaryList$halfGenomeRange$zoneCol_alpha,
		inherit.aes = FALSE) +
	guides(colour = guide_legend(
		title = "Gene Group",
		label.position = "right",
		title.position = "top")) +
	coord_cartesian(xlim = c(0, 0.5), expand = FALSE) +
	lightTheme + 
	theme(
		panel.spacing = unit(2, "lines"),
		panel.grid.minor = element_blank(),
		legend.position = "right",
		legend.background = element_rect(color = axisCol, fill = "white"),
		legend.text.align = 0.5,
		legend.text = element_text(margin = margin(t = 5, b = 5)),
		strip.background = element_rect(color = axisCol, fill = "white"),
		strip.text = element_text(color = textCol, size = 12)
	)

quartz(width = 18, height = 8)
print(firstThird_orth_plot)
quartz.save(file = file.path(orthPositionFig_extPath, paste0(testCOG, "_All_firstThird_density.pdf")), type = "pdf", dpi = 300)
invisible(dev.off())

# -------------------------------------- #
# Other two thirds
twoThird_genes	<- GPA_testCOG_byType %>% filter(distToOri > nadirDist & type == "All")


twoThird_orthologs <- processOrthologPosition_2(
	COGname = testCOG,
	orthData = twoThird_genes,
	dbConn = dbConn,
	dnaA_data = dnaA_clean_trim,
	speciesSubgroups = speciesSubgroups_df,
	superGroupTaxid_list = superGroup_list,
	outlierTaxid = outlierTaxid)


twoThird_posBias_dt	<- as.data.table(twoThird_orthologs$supergroupCompareData)
twoThird_posBias_dt	<- twoThird_posBias_dt %>% 
	mutate(sporeForming = case_when(
			SuperGroup %in% sporeGroups ~ "Spore",
			SuperGroup %in% nonSporeGroups ~ "No Spore",
			TRUE ~ "NA")
	) %>% 
	mutate(sporeForming = factor(sporeForming, levels = c("Spore", "No Spore")))

# Select those entries where Orths are found in both Spore and Non-spore and those where only in one
twoThird_posBias_dt %<>%
	group_by(Group) %>% 
	mutate(sporePresence = case_when(
		n_distinct(sporeForming) == 2 ~ "Both",
		sporeForming == "Spore" ~ "S Only",
		TRUE ~ "NS Only"
		)
	)

# Extend the dataset for plotting
twoThird_posBias_dt_ext	<-  twoThird_posBias_dt %>%
	subset(distToOri < 0.1 | distToOri > 0.4) %>%
	mutate(distToOri = case_when(
		distToOri < 0.1 ~ distToOri + 0.5,
		TRUE ~ distToOri - 0.5)
	) %>%
	bind_rows(twoThird_posBias_dt)


# Plot density of genes in groups that have Orths restricted to either S or non-S only
speciesRestrcited_twoThird <- twoThird_posBias_dt_ext %>% group_by(Group, sporePresence) %>% filter(sporePresence != "Both")
speciesRestrcited_twoThird_plot <- ggplot(data = speciesRestrcited_twoThird, mapping = aes(x = distToOri, y = ..count.., color = sporePresence, fill = sporePresence)) +
	scale_y_continuous(name = "Gene density") +
	scale_x_continuous(name = "Relative distance from origin") +
	stat_density(geom = "line", adjust = 1/10, size = 2, position = "identity") +
	scale_color_manual(values = sporeCols) +
	scale_fill_manual(values = alpha(sporeCols, 0), guide = FALSE) +
	geom_rect(
		data = zoneBoundaryList$halfGenomeRange,
		aes(xmin = zoneMin, xmax = zoneMax, ymin = -Inf, ymax = Inf),
		fill = zoneBoundaryList$halfGenomeRange$zoneCol_alpha,
		inherit.aes = FALSE) +
	guides(colour = guide_legend(
		title = "Gene Group",
		label.position = "right",
		title.position = "top")) +
	coord_cartesian(xlim = c(0, 0.5), expand = FALSE) +
	lightTheme + 
	theme(
		panel.spacing = unit(2, "lines"),
		panel.grid.minor = element_blank(),
		legend.position = "right",
		legend.background = element_rect(color = axisCol, fill = "white"),
		legend.text.align = 0.5,
		legend.text = element_text(margin = margin(t = 5, b = 5)),
		strip.background = element_rect(color = axisCol, fill = "white"),
		strip.text = element_text(color = textCol, size = 12)
	)

# quartz(width = 18, height = 8)
# print(speciesRestrcited_twoThird_plot)
# quartz.save(file = file.path(orthPositionFig_extPath, paste0(testCOG, "_TwoThird_speciesRestricted_density.pdf")), type = "pdf", dpi = 300)
# invisible(dev.off())


twoThird_orth_plot	<- ggplot(data = twoThird_posBias_dt_ext, mapping = aes(x = distToOri, y = ..count.., color = sporeForming, fill = sporeForming)) +
	scale_y_continuous(name = "Gene density") +
	scale_x_continuous(name = "Relative distance from origin") +
	stat_density(geom = "line", adjust = 1/10, size = 2, position = "identity") +
	scale_color_manual(values = sporeCols) +
	scale_fill_manual(values = alpha(sporeCols, 0), guide = FALSE) +
	geom_rect(
		data = zoneBoundaryList$halfGenomeRange,
		aes(xmin = zoneMin, xmax = zoneMax, ymin = -Inf, ymax = Inf),
		fill = zoneBoundaryList$halfGenomeRange$zoneCol_alpha,
		inherit.aes = FALSE) +
	guides(colour = guide_legend(
		title = "Gene Group",
		label.position = "right",
		title.position = "top")) +
	coord_cartesian(xlim = c(0, 0.5), expand = FALSE) +
	lightTheme + 
	theme(
		panel.spacing = unit(2, "lines"),
		panel.grid.minor = element_blank(),
		legend.position = "right",
		legend.background = element_rect(color = axisCol, fill = "white"),
		legend.text.align = 0.5,
		legend.text = element_text(margin = margin(t = 5, b = 5)),
		strip.background = element_rect(color = axisCol, fill = "white"),
		strip.text = element_text(color = textCol, size = 12)
	)
quartz(width = 18, height = 8)
print(twoThird_orth_plot)
quartz.save(file = file.path(orthPositionFig_extPath, paste0(testCOG, "_All_twoThird_density.pdf")), type = "pdf", dpi = 300)
invisible(dev.off())











