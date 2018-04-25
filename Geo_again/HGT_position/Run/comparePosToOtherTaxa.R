#!/usr/bin/env Rscript

# Load master variables and HGT position functions
invisible(sapply(HGTPos.all, source, .GlobalEnv))

# Any libraries explicitly used in the script
require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load("RSQLite", "dplyr","ggplot2", "ggpubr", "wesanderson")


# ------------------------------------------------------------------------------------- #
# Read in data

message("\nReading in data...", appendLF = FALSE)

# General position data
perTypeData			<- readRDS(file.path(positionData_path, "AG_perTypeData.rds"))

# Per-type COG data
perTypeCOG_data		<- readRDS(file.path(positionData_path, "AG_perTypeCOGData.rds"))

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
byCOG_geneSpreadSummary_df$COG	<- factor(byCOG_geneSpreadSummary_df$COG, levels = levels(eventNumberByCOG$COG))

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
	event_B		<- sample(unique(subset(byCOG_geneSpreadAll_df, COGcat %in% event_COG, select = eventIndex, drop = TRUE)), 1)

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
# Bind the within-Group gene coLocation with the random sampling
compToRandom_df	<- bind_rows(list(subset(perGroupAlt_df, select = c(range, type)), rand_df, randInCOG_df))
compToRandom_df$type	<- factor(compToRandom_df$type, levels = c("Random", "RandomInCOG", "IntraGroup"))

# We'll do a wilcox test on all combinations (3)
statComparisons		<- lapply(combn(unique(compToRandom_df$type), 2, simplify = FALSE), paste0)
plotCols			<- wes_palette("Rushmore1")[3:5]

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
speciesSubgroups_df		<- read.table(speciesSubgroups_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

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




hgtCOGCat_list	<- unlist(unique(soloEventData_list$COGcat))

lapply(hgtCOGCat_list, function(COGname) {

	onlyCOG_data			<- subset(soloEventData_list, COGcat == COGname)
	onlyCOG_data$dist2ori	<- ifelse(onlyCOG_data$relGeneStart > 0.5, 1 - onlyCOG_data$relGeneStart, onlyCOG_data$relGeneStart)
	onlyCOG_groups			<- unique(onlyCOG_data$orthGroup)



	byGroupOrthPosition	<- lapply(onlyCOG_groups, function(group) {

		message(paste0("COG == ", COGname, ": working on group ", group))

		## All the data corresponding to the orthologs
		orthGroup_tbl	<- dbSendQuery(dbConn, 'SELECT * FROM t1 WHERE orthGroup = :orthGroup')
		dbBind(orthGroup_tbl, param = list(orthGroup = group))
		orthGroup_df	<- dbFetch(orthGroup_tbl)
		dbClearResult(orthGroup_tbl)

		## Select only those that have matching taxids in the clean DNAa_data
		pruneOrth_data	<- orthGroup_df[which(orthGroup_df$taxid %in% dnaA_clean_df$taxid),]
		# Select only those that are on the main chromosome
		pruneOrth_data	<- pruneOrth_data[which(pruneOrth_data$plasmid == "F"),]
		# Clean memory
		rm(orthGroup_df)

		# ------------------------------------------------------------------------------------- #

	
		combineOri_df	<- left_join(pruneOrth_data, dnaA_clean_trim, by = "taxid")
		combineOri_df$relGeneStart	<- apply(combineOri_df, 1, genomeRelativePosition_format)
		combineOri_df$dist2ori		<- ifelse(combineOri_df$relGeneStart > 0.5, 1 - combineOri_df$relGeneStart, combineOri_df$relGeneStart)
		combineOri_df	<- subset(combineOri_df, select = -c(sequence, NuclSeq))

		perGroupGenPos_df <- combineOri_df
	
		# WIP
		# WIP
		# WIP
		# WIP
		# WIP
		# WIP
		# WIP
		# ------------------------------------------------------------------------------------- #

		## Summarise by subgroup
		clustB			<- makeCluster(numCores, type = "FORK")
		bySubgroup_df	<- summariseSubgroups(df = perGroupGenPos_df, subgroupData = speciesSubgroups_df, summaryStat = "dist2ori", clusterCon = clustB)
		stopCluster(clustB)
		
		# ------------------------------------------------------------------------------------- #

		## Gene data for the HGT only for this group
		HGT_entries		<- perTypeData$lHGT$`4`$allPosData$protID[which(perTypeData$lHGT$`4`$allPosData$orthGroup == group)]
		HGT_only_set	<- perGroupGenPos_df[which(perGroupGenPos_df$protID %in% HGT_entries),]

		plotTemp	<- ggplot(data = perGroupGenPos_df, aes(dist2ori)) + 
			geom_density(adjust = 1/5) +
			geom_density(data = bySubgroup_df, adjust = 1/5, color = "red") +
			geom_histogram(
				data = subset(perGroupGenPos_df, is_ag == 1),
				aes(x = dist2ori, y = ..ncount..),
				bins = 100, 
				fill = "red",
				inherit.aes = FALSE) +
			geom_histogram(
				data = HGT_only_set,
				aes(x = dist2ori, y = -..ncount..),
				bins = 100, 
				fill = "blue",
				inherit.aes = FALSE)

		return(list(allData = bySubgroup_df, AG_HGT_data = HGT_only_set, numberOfGenes = nrow(perGroupGenPos_df), numberOfSpecies = nrow(bySubgroup_df), plot = plotTemp))
	})

	## Rename the list by group
	names(byGroupOrthPosition)	<- paste0(COGname, onlyCOG_groups)

})








byGroupOrthPosition_df	<- bind_rows(lapply(byGroupOrthPosition, function(element) return(element$allData)))
noPlasmid				<- byGroupOrthPosition_df[which(byGroupOrthPosition_df$plasmid == "F"),]


bacillaceae_accAss	<- read.table(file = file.path(master_dir, "Consensus_groups", "Bacillaceae", "Bacillaceae_acc_ass_list.txt"), sep = "\n", header = FALSE, stringsAsFactors = FALSE)$V1
bacillales_accAss	<- read.table(file = file.path(master_dir, "Consensus_groups", "Bacillales", "Bacillales_acc_ass_list.txt"), sep = "\n", header = FALSE, stringsAsFactors = FALSE)$V1
clostridia_accAss	<- read.table(file = file.path(master_dir, "Consensus_groups", "Clostridia", "Clostridia_acc_ass_list.txt"), sep = "\n", header = FALSE, stringsAsFactors = FALSE)$V1

AGOnly			<- noPlasmid[which(noPlasmid$binomial %in% binomial_list),]
bacillalesOnly	<- noPlasmid[which(noPlasmid$acc_ass %in% bacillales_accAss),]
bacillaceaeOnly	<- noPlasmid[which(noPlasmid$acc_ass %in% bacillaceae_accAss),]
clostridiaOnly	<- noPlasmid[which(noPlasmid$acc_ass %in% clostridia_accAss),]

Orth_G_AG_plot			<- ggplot(data = AGOnly, aes(dist2ori)) + geom_density(adjust = 1/4, col = "red") + geom_density(data = justG, adjust = 1/4, col = "blue")
Orth_G_Bacillaceae_plot	<- ggplot(data = bacillaceaeOnly, aes(dist2ori)) + geom_density(adjust = 1/4)
Orth_G_Bacillales_plot	<- ggplot(data = bacillalesOnly, aes(dist2ori)) + geom_density(adjust = 1/4)
Orth_G_Clostridia_plot	<- ggplot(data = clostridiaOnly, aes(dist2ori)) + geom_density(adjust = 1/4)







dbDisconnect(conn)

