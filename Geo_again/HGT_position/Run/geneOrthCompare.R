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

# Solo HGT event dataframe
soloEventData_df	<- readRDS(file.path(positionData_path, "soloHGTEvent_data.rds"))#

# Zone list
zoneBoundaryList	<- readRDS(file.path(positionData_path, "AG_zoneBoundaries.rds"))

message("\rReading in data... done\n")

# -------------------------------------- #
# Output path for figures
orthPositionFig_path	<- file.path(figureOutput_path, "orthPosition")
if(!dir.exists(orthPositionFig_path)) dir.create(orthPositionFig_path)



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

sporeGroups		<- c("Bacillaceae", "Paenibacillaceae", "Alicyclobacillaceae", "Pasteuriaceae", "Sporolactobacillaceae")
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

# -------------------------------------- #

sporeCols	<- c(wes_palette("Chevalier1")[1:2], wes_palette("FantasticFox1")[5])

cutOffCogs	<- all_posBias_dt %>% group_by(COGcat) %>% filter(sporeForming == "HGT") %>% summarise(count = n()) %>% filter(count >= 100) %>% pull(COGcat)
all_posBias_numCut_dt	<- all_posBias_dt %>% filter(COGcat %in% cutOffCogs)

allCOGs_sporeCompare_plot	<- ggplot(data = all_posBias_numCut_dt, mapping = aes(x = distToOri, color = sporeForming, fill = sporeForming)) +
	facet_wrap(~COGcat, scales = "free_y") +
	scale_y_continuous(name = "Gene density") +
	scale_x_continuous(name = "Relative distance from origin") +
	stat_density(geom = "line", adjust = 1/8, size = 0.7, position = "identity") +
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



quartz(width = 24, height = 12)
print(allCOGs_sporeCompare_plot)
quartz.save(file = file.path(orthPositionFig_path, "allCOG_sporeVSnonspore.pdf"), type = "pdf", dpi = 300)
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



