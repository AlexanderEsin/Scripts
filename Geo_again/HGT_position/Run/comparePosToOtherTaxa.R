#!/usr/bin/env Rscript

# Load master variables and HGT position functions
invisible(sapply(HGTPos.all, source, .GlobalEnv))

# Any libraries explicitly used in the script
require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load("RSQLite", "dplyr","ggplot2", "ggpubr", "wesanderson", "gridExtra")


# ------------------------------------------------------------------------------------- #
# Read in data

message("\nReading in data...", appendLF = FALSE)

# General position data
perTypeData			<- readRDS(file.path(positionData_path, "AG_perTypeData.rds"))

# Per-type COG data
perTypeCOG_data		<- readRDS(file.path(positionData_path, "AG_perTypeCOGData.rds"))

# taxdmp nodes & names data
nodes_data			<- readRDS(file.path(taxdmp_path, "nodes.rds"))
names_data			<- readRDS(file.path(taxdmp_path, "names.rds"))

# Subdivision key
subDivisionKey_file	<- file.path(positionData_path, "subDivisionKeyData.rds")
if(!file.exists(subDivisionKey_file)) {
	stop("The subdivision data file \"", subDivisionKey_file, "\" is required.\nRun the subDivisionPrep.R script!")
} else {
	subDivisionKey_data	<- readRDS(file.path(positionData_path, "subDivisionKeyData.rds"))
}

message("\rReading in data... done\n")

# ------------------------------------------------------------------------------------- #
# Open All_prot database
dbConn			<- dbConnect(RSQLite::SQLite(), allProtDB_path)


# ------------------------------------------------------------------------------------- #

# List of all the unique COG categorties found in the entire AG dataset
uniqueCOGs		<- unique(unlist(perTypeData$All$allPosData$COGcat))

# For each COG and each individual HGT event belonging to that COG, find the range of distance to origin
# Calculate the max range and SD of the gene positions per-event/per-COG
byCOG_geneSpread_list	<- lapply(uniqueCOGs, function(COG) {

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
		eventDistSD		<- sd(oriDists)
		eventDistMax	<- max(oriDists) - min(oriDists)

		# All members of event must belong to single Orth group. Sanity check
		orthGroup	<- unique(perEvent_data$orthGroup)
		if (length(orthGroup) > 1) {
			stop("Should not have more than one group for HGT event")
		}

		# Drop the CircStart and CircEnd columns, because combining circular vectors with bind_rows throws warnings
		perEvent_data	<- subset(perEvent_data, select = -c(CircStart, CircEnd))
		perEvent_data$dist2ori	<- oriDists

		# Summary df
		summary_df	<- data.frame(COG = COG, orthGroup = orthGroup, eventIndex = hgtEvent, totalGenes = numTotal, geneSD = eventDistSD, geneDistMax = eventDistMax, stringsAsFactors = FALSE)
		return(list(allData = perEvent_data, summary_df = summary_df))
	})

	# Combine the per-event results (into per-COG results)
	byEventAllData_df 	<- bind_rows(lapply(byEventSpread_list, function(element) return(element$allData)))
	byEventSummary_df 	<- bind_rows(lapply(byEventSpread_list, function(element) return(element$summary_df)))
	return(list(allData = byEventAllData_df, summaryData = byEventSummary_df))
})
# Combine the cross-COG results
byCOG_geneSpreadAll_df		<- bind_rows(lapply(byCOG_geneSpread_list, function(element) return(element$allData)))
byCOG_geneSpreadSummary_df	<- bind_rows(lapply(byCOG_geneSpread_list, function(element) return(element$summaryData)))


# Which COGs have the most independent events?
eventNumberByCOG	<- sort(table(byCOG_geneSpreadSummary_df$COG), decreasing = TRUE)
eventNumberByCOG_df	<- data.frame(eventNumberByCOG)
names(eventNumberByCOG_df)	<- c("COG", "Freq")

# Sort the COG categories by the number of events
byCOG_geneSpreadSummary_df$COG	<- factor(byCOG_geneSpreadSummary_df$COG, levels = levels(eventNumberByCOG_df$COG))

# Plot the maximum distance between genes (relative to distance from Ori) for every lHGT event
byMaxGeneDist_perEvent_boxplot	<- ggplot(data = subset(byCOG_geneSpreadSummary_df, totalGenes > 1), aes(x = COG, y = geneDistMax, group = COG)) +
	geom_boxplot(fill = "grey40", color = "#D9D9D9") +
	geom_label(data = eventNumberByCOG_df, aes(x = COG, y = 0.5, label = Freq), fill = "#333233", color = "#D9D9D9", inherit.aes = FALSE) +
	darkTheme

# Plot the SD of the distance between genes (relative to distance from Ori) for every lHGT event
byGeneSD_perEvent_boxplot	<- ggplot(data = subset(byCOG_geneSpreadSummary_df, totalGenes > 1), aes(x = COG, y = geneSD, group = COG)) +
	geom_boxplot(fill = "grey40", color = "#D9D9D9") +
	geom_label(data = eventNumberByCOG_df, aes(x = COG, y = 0.3, label = Freq), fill = "#333233", color = "#D9D9D9", inherit.aes = FALSE) +
	darkTheme



# ------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------- #
# Check which orthologous groups occur more than once - indicating multiple events in same gene family
multiEventGroups_list	<- unique(byCOG_geneSpreadSummary_df$orthGroup[duplicated(byCOG_geneSpreadSummary_df$orthGroup)])

# For each gene family that contains more than one HGT event, find the mean position per event
# Compare this to the other events in the same gene family.
# Naive assumption: independent events in same gene family should be positioned nearby
perGroup_list	<- lapply(multiEventGroups_list, function(orthGroup) {
	perGroup_data	<- byCOG_geneSpreadAll_df[which(byCOG_geneSpreadAll_df$orthGroup == orthGroup),]
	perGroup_events	<- unique(perGroup_data$eventIndex)

	perEvent_list	<- lapply(perGroup_events, function(hgtEvent) {
		perEvent_data	<- perGroup_data[which(perGroup_data$eventIndex == hgtEvent),]
		oriDists		<- unlist(lapply(perEvent_data$relGeneStart, function(geneStart) ifelse(geneStart > 0.5, 1 - geneStart, geneStart)))
		perEvent_sd		<- sd(oriDists)
		perEvent_mean	<- mean(oriDists)
		out_df			<- data.frame(orthGroup = orthGroup, eventIndex = hgtEvent, COGcat = unique(unlist(perEvent_data$COGcat)), numGenes = nrow(perEvent_data), meanLocation = perEvent_mean, sdLocation = perEvent_sd, stringsAsFactors = FALSE)
		return(out_df)
	})
	perEvent_combined	<- bind_rows(perEvent_list)
	perEvent_combined$range	<- max(perEvent_combined$meanLocation) - min(perEvent_combined$meanLocation)
	return(perEvent_combined)
})

# Extract the group Number, COGcat, and the maximum distance between the mean per-Event positions
perGroupAlt		<- lapply(perGroup_list, subset, c(TRUE, rep(FALSE, 100)), select = c(orthGroup, COGcat, range))
perGroupAlt_df	<- bind_rows(perGroupAlt)

# Order the df ascending in range
perGroupAlt_df	<- perGroupAlt_df[order(perGroupAlt_df$range),]
perGroupAlt_df$plotIndex	<- seq(1:nrow(perGroupAlt_df))

# Which COGs have the most independent events?
groupNumberByCOG	<- sort(table(perGroupAlt_df$COGcat), decreasing = TRUE)
groupNumberByCOG	<- data.frame(groupNumberByCOG)
names(groupNumberByCOG)	<- c("COGcat", "Freq")

# Sort the COG categories by the number of events
perGroupAlt_df$COGcat	<- factor(perGroupAlt_df$COGcat, levels = levels(groupNumberByCOG$COGcat))
perGroupAlt_df$type		<- "IntraGroup"

# ------------------------------------------------------------------------------------- #
# Random sampling to see whether independent HGT events within the same orthologous group co-located more or less than expected

numObserv	<- nrow(perGroupAlt_df)

# Random comparison 1: the distance between the mean positions of genes in any HGT event A and any HGT event B
events_A	<- sample(unique(byCOG_geneSpreadAll_df$eventIndex), numObserv)
events_B	<- sample(unique(byCOG_geneSpreadAll_df$eventIndex), numObserv)
randEventRange	<- unlist(lapply(1:length(events_A), function(eventNum) {
	eventA	<- events_A[eventNum]
	eventB	<- events_B[eventNum]

	eventComparison	<- unlist(lapply(c(eventA, eventB), function(event) {
		event_data		<- byCOG_geneSpreadAll_df[which(byCOG_geneSpreadAll_df$eventIndex == event),]
		event_oriDist	<- unlist(lapply(event_data$relGeneStart, function(geneStart) ifelse(geneStart > 0.5, 1 - geneStart, geneStart)))
		event_mean		<- mean(event_oriDist)
		return(event_mean)
	}))

	xEventRange		<- max(eventComparison) - min(eventComparison)
	return(xEventRange)
}))
rand_df	<- data.frame(range = randEventRange, type = "Random", stringsAsFactors = FALSE)


# Random comparison 2: as above, but event B must belong to the same COGcat as event A
events_A	<- sample(unique(byCOG_geneSpreadAll_df$eventIndex), numObserv)
randInCOGEventRange	<- unlist(lapply(1:length(events_A), function(eventNum) {
	eventA		<- events_A[eventNum]
	event_COG	<- unique(unlist(byCOG_geneSpreadAll_df$COG[which(byCOG_geneSpreadAll_df$eventIndex == eventA)]))
	print(event_COG)
	eventB		<- sample(unique(subset(byCOG_geneSpreadAll_df, COGcat == event_COG, select = eventIndex, drop = TRUE)), 1)

	eventComparison	<- unlist(lapply(c(eventA, eventB), function(event) {
		event_data		<- byCOG_geneSpreadAll_df[which(byCOG_geneSpreadAll_df$eventIndex == event),]
		event_oriDist	<- unlist(lapply(event_data$relGeneStart, function(geneStart) ifelse(geneStart > 0.5, 1 - geneStart, geneStart)))
		event_mean		<- mean(event_oriDist)
		return(event_mean)
	}))

	xEventRange		<- max(eventComparison) - min(eventComparison)
	return(xEventRange)
}))
randInCOG_df	<- data.frame(range = randInCOGEventRange, type = "RandomInCOG", stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------------- #
# Partition the multiEvent data further by the COG they represent - some COGs have an intrinsic positional bias

# Positionally biased and non-biased COGs
positionBiasCOGs		<- c("J", "M", "G", "C", "E", "Q", "O")
positionNonBiasCOGs		<- c("T", "K", "V", "L", "P", "I")

# Partition and label data
perGroupBiasCOGs		<- perGroupAlt_df[which(perGroupAlt_df$COGcat %in% positionBiasCOGs),]
perGroupNonBias			<- perGroupAlt_df[which(perGroupAlt_df$COGcat %in% positionNonBiasCOGs),]
perGroupBiasCOGs$type	<- "IntraCOGBias"
perGroupNonBias$type	<- "IntraCOGNonBias"

# Bind the random sample data with multi event data
compToRandom_df			<- bind_rows(list(
	subset(perGroupAlt_df, select = c(range, type)),
	subset(perGroupBiasCOGs, select = c(range, type)),
	subset(perGroupNonBias, select = c(range, type)),
	rand_df,
	randInCOG_df))

# Factor the type for plot order
typeOrder				<- c("Random", "RandomInCOG", "IntraGroup", "IntraCOGBias", "IntraCOGNonBias")
compToRandom_df$type	<- factor(compToRandom_df$type, levels = typeOrder)

# We'll do a wilcox test on all combinations of random and full multiEvent set + 1 comparision between positional and non-positioned COGs
statComparisons_1	<- lapply(combn(typeOrder[1:3], 2, simplify = FALSE), paste0)
statComparisons_2	<- lapply(combn(typeOrder[4:5], 2, simplify = FALSE), paste0)
statComparisons		<- c(statComparisons_1, statComparisons_2)

# Plot colours
plotCols			<- wes_palette("Rushmore1")[1:5]

# Produce the boxplot
multiInGroupEvents_boxplot <- ggplot(data = compToRandom_df, aes(x = type, y = range, fill = type)) +
	scale_y_continuous(name = "Max range in mean location of HGT events") +
	geom_boxplot(color = "#D9D9D9") +
	scale_fill_manual(values = plotCols, guide = FALSE) +
	stat_compare_means(mapping = aes(x = type, y = range), comparisons = statComparisons, method = "wilcox.test", p.adjust = "bonferroni", color = "#D9D9D9", size = 0.5, label = "p.format") +
	darkTheme +
	theme(
		axis.title.x = element_blank(),
		panel.grid.major.x = element_blank()
	)



# ------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------- #

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


# ------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------- #

## Entries where there is only a single HGT event per orthGroup
soloEventGroups_list	<- byCOG_geneSpreadSummary_df[which(!byCOG_geneSpreadSummary_df$orthGroup %in% multiEventGroups_list),]
soloEventData_list		<- byCOG_geneSpreadAll_df[which(byCOG_geneSpreadAll_df$orthGroup %in% soloEventGroups_list$orthGroup),]


byCOGDispersionSolo_l	<- lapply(unique(unlist(soloEventData_list$COGcat)), function(COG) {

	byCOGsoloEvent_data	<- soloEventData_list[which(soloEventData_list$COGcat == COG),]
	uniqueGroups		<- unique(byCOGsoloEvent_data$orthGroup)

	perCOGbyGroup_list	<- lapply(uniqueGroups, function(group) {
		byGroup_data	<- byCOGsoloEvent_data[which(byCOGsoloEvent_data$orthGroup == group),]
		oriDists		<- unlist(lapply(byGroup_data$relGeneStart, function(geneStart) ifelse(geneStart > 0.5, 1 - geneStart, geneStart)))
		meanLocation	<- mean(oriDists)
		sdLocation		<- sd(oriDists)
		out_df			<- data.frame(orthGroup = group, COGcat = COG, numGenes = nrow(byGroup_data), meanLocation = meanLocation, sdLocation = sdLocation, stringsAsFactors = FALSE)
		return(out_df)
	})

	perCOGbyGroup_df	<- bind_rows(perCOGbyGroup_list)
	return(perCOGbyGroup_df)
})
byCOGDispersionSolo_df	<- bind_rows(byCOGDispersionSolo_l)

# ------------------------------------------------------------------------------------- #
# SuperGroup determination - we want to be able to filter orthologs by some taxonomic groups
superGroups		<- c("Bacillaceae", "Clostridia", "Bacillales")

recalculateSGroups		<- TRUE
rdsSuperGroup_file		<- file.path(positionData_path, "parentSpecies_bySupergroup.rds")
if (file.exists(rdsSuperGroup_file)) {
	bySupergroupTaxid_list	<- readRDS(rdsSuperGroup_file)
	taxaInFile	<- names(bySupergroupTaxid_list)
	
	if (identical(superGroups[order(superGroups)], taxaInFile[order(taxaInFile)])) {
		message("Taxonomic groupings successfully read from file")
		recalculateSGroups		<- FALSE
	}
}

if (recalculateSGroups) {
	message("Recalculating taxonomic groupings")
	# Identify the parentSpecies taxids that belong to desired superGroups
	bySupergroupTaxid_list	<- lapply(superGroups, taxidsBySupergroup, speciesGroupings = speciesSubgroups_df, nodes_data = nodes_data, names_data = names_data)
	names(bySupergroupTaxid_list)	<- superGroups
	rm(list = c("names_data", "nodes_data"))

	# Save object
	saveRDS(object = bySupergroupTaxid_list, file = file.path(positionData_path, "parentSpecies_bySupergroup.rds"))
}

# ------------------------------------------------------------------------------------- #

# All COG categories in all solo event Data
allSoloCOGCat	<- unlist(unique(soloEventData_list$COGcat))

# All ortholog data
perCOGOrth_data			<- lapply(allSoloCOGCat, processOrthologPosition,
	orthData = soloEventData_list,
	dbConn = dbConn,
	superGroupTaxid_list = bySupergroupTaxid_list)
names(perCOGOrth_data)	<- allSoloCOGCat


# Orthologs of AG HGT groups, of which all members are in Ori/Near-Ori
perCOGOriOrth_data		<- lapply(allSoloCOGCat, processOrthologPosition,
	orthData = soloEventData_list,
	distRange = c(0, 0.25),
	dbConn = dbConn,
	superGroupTaxid_list = bySupergroupTaxid_list
)
names(perCOGOriOrth_data)	<- nearOriCOGCat

# Orthologs of AG HGT groups, of which all members are in Ter/Near-Ter
perCOGTerOrth_data		<- lapply(allSoloCOGCat, processOrthologPosition,
	orthData = soloEventData_list,
	distRange = c(0.25, 0.5),
	dbConn = dbConn,
	superGroupTaxid_list = bySupergroupTaxid_list
)
names(perCOGTerOrth_data)	<- nearOriCOGCat


# ------------------------------------------------------------------------------------- #
# Define the ORI-enriched and TER-enriched functional COGs
oriBiasCOGs		<- c("J", "M")
terBiasCOGs		<- c("C", "E", "G", "Q")

# Overall co-occurence of ortholog genes in the same locations?
all_terBias_data	<- bind_rows(lapply(terBiasCOGs, function(COG) return(perCOGOrth_data[[COG]]$supergroupCompareData)))
all_oriBias_data	<- bind_rows(lapply(oriBiasCOGs, function(COG) return(perCOGOrth_data[[COG]]$supergroupCompareData)))
all_posBias_data	<- bind_rows(list(all_oriBias_data, all_terBias_data))
all_posBias_data$COGcat	<- factor(all_posBias_data$COGcat, levels = c(oriBiasCOGs, terBiasCOGs))

# Do metabolic genes that are usually near Ori or near Ter have co-positioned orthologs?
ter_terBias_data	<- bind_rows(lapply(terBiasCOGs, function(COG) return(perCOGTerOrth_data[[COG]]$supergroupCompareData)))
ori_terBias_data	<- bind_rows(lapply(terBiasCOGs, function(COG) return(perCOGOriOrth_data[[COG]]$supergroupCompareData)))

# Set up colors
allCol	<- wes_palette("Chevalier1")[3]
grpCol	<- wes_palette("FantasticFox1")[1:length(superGroups)]
hgtCol	<- wes_palette("Darjeeling1")[1]

# ------------------------------------------------------------------------------------- #
# Overall co-occurence of ortholog genes
positionBiasOrthCompare_plot	<- ggplot(data = all_posBias_data, mapping = aes(x = dist2ori, color = SuperGroup)) +
	facet_wrap(~COGcat, ncol = 2, nrow = 3, scales = "free_y") +
	scale_y_continuous(name = "Gene density") +
	scale_x_continuous(name = "Relative distance from origin") +
	stat_density(geom = "line", adjust = 1/5, position = "identity") +
	scale_color_manual(values = c(grpCol, hgtCol), name = "Gene Group") +
	darkTheme + 
	theme(
		panel.grid.minor.x = element_blank(),
		legend.background = element_rect(color = "#D9D9D9", fill = "#333233"),
		strip.background = element_rect(color = "#D9D9D9", fill = "#333233"),
		strip.text = element_text(color = "#D9D9D9", size = 12)
	)

Ter_MetabolicOrth_plot		<- ggplot(data = ter_terBias_data, mapping = aes(x = dist2ori, color = SuperGroup)) +
	facet_wrap(~COGcat, nrow = 2, scales = "free") +
	scale_y_continuous(name = "Gene density") +
	scale_x_continuous(name = "Relative distance from origin") +
	stat_density(geom = "line", adjust = 1/5, position = "identity") +
	scale_color_manual(values = c(grpCol, hgtCol), name = "Gene Group") +
	darkTheme + 
	theme(
		legend.background = element_rect(color = "#D9D9D9", fill = "#333233"),
		strip.background = element_rect(color = "#D9D9D9", fill = "#333233"),
		strip.text = element_text(color = "#D9D9D9", size = 12)
	)

Ori_MetabolicOrth_plot	<- ggplot(data = ori_terBias_data, mapping = aes(x = dist2ori, color = SuperGroup)) +
	facet_wrap(~COGcat, nrow = 2, scales = "free") +
	scale_y_continuous(name = "Gene density") +
	scale_x_continuous(name = "Relative distance from origin") +
	stat_density(geom = "line", adjust = 1/5, position = "identity") +
	scale_color_manual(values = c(grpCol, hgtCol), name = "Gene Group") +
	darkTheme + 
	theme(
		legend.background = element_rect(color = "#D9D9D9", fill = "#333233"),
		strip.background = element_rect(color = "#D9D9D9", fill = "#333233"),
		strip.text = element_text(color = "#D9D9D9", size = 12)
	)





dbDisconnect(conn)

