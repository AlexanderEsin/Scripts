#!/usr/bin/env Rscript

# Load master variables and HGT position functions
invisible(sapply(HGTPos.all, source, .GlobalEnv))

# Any libraries explicitly used in the script
require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load("RSQLite", "dplyr","ggplot2", "ggpubr", "wesanderson", "gridExtra", "stringr")


# ------------------------------------------------------------------------------------- #
message("\nReading in data...", appendLF = FALSE)

# General position data
perTypeData			<- readRDS(file.path(positionData_path, "AG_perTypeData.rds"))

message("\rReading in data... done\n")

# Open All_prot database
dbConn			<- dbConnect(RSQLite::SQLite(), allProtDB_path)

# ------------------------------------------------------------------------------------- #
# Extend the data an extra 25% to cover the origin fully (remove edge effects in density estimation)
dataTypes_trunc		<- c("All", "lHGT", "Ver")
withExtendedPos_data <- lapply(dataTypes_trunc, function(dataType) {

	# Identify the dataframe we need for the dataType
	if (identical(dataType, "All")) {
		data	<- 	perTypeData$All$allPosData
	} else if (identical(dataType, "Ver")) {
		data	<- perTypeData[[dataType]][[verPenalty]]$allPosData
	} else {
		data	<- perTypeData[[dataType]][[hgtPenalty]]$allPosData
	}

	# Subset to remove the circular start and end
	data_trim	<- subset(data, select = -c(CircStart, CircEnd))
	
	# Add the dataType as a column
	data_trim$type	<- dataType

	# Add the first quarter of the data onto the end to calculate densities over the Origin
	firstQ		<- subset(data_trim, relGeneStart >= 0 & relGeneStart <= 0.25)
	firstQ$relGeneStart	<- firstQ$relGeneStart + 1
	firstQ$relGeneEnd	<- firstQ$relGeneEnd + 1
	extended	<- bind_rows(list(data_trim, firstQ))

	return(list(normal_df = data_trim, extended_df = extended))
})

# Without extension
normalPos_df		<- bind_rows(lapply(withExtendedPos_data, function(x) return(x$normal_df)))
normalPos_df$type	<- factor(normalPos_df$type, levels = dataTypes_trunc)

# With 25% extension of data
extendedPos_df		<- bind_rows(lapply(withExtendedPos_data, function(x) return(x$extended_df)))
extendedPos_df$type	<- factor(extendedPos_df$type, levels = dataTypes_trunc)

# ------------------------------------------------------------------------------------- #
# Exploratory plot of linear HGT vs Ver densities to identify boundaries

# Set density curve colours
all_HGT_Ver_cols	<- c(alpha(wes_palette("Chevalier1")[4], 0.5), wes_palette("Darjeeling1")[2:3])

# Plot
linearXSpecies_verHGT_dens_plot <- ggplot(data = extendedPos_df, aes(x = relGeneStart, color = type)) +
	scale_x_continuous(
		name = "Normalized genome position",
		limits = c(0, 1.25),
		breaks = seq(0, 1, by = 0.5),
		labels = c("Origin", "Terminus", "Origin"),
		minor_breaks = seq(0, 1.25, by = ((1 / 360) * 30)),
		sec.axis = dup_axis(
			name = NULL,
			breaks = seq(0, 1.25, by = ((1 / 360) * 30)),
			labels = c("Origin", seq(30, 150, by = 30), "Terminus", seq(210, 330, by = 30) ,"Origin", seq(30, 90, by = 30))
		)) +
	scale_y_continuous(name = "Gene Enrichment") +
	stat_density(geom = "line", position = "identity", n = 2^12, adjust = 1/10, size = 1) +
	scale_color_manual(values = all_HGT_Ver_cols) +
	ggtitle("Exploratory zone boundary analysis") +
	theme_classic() +
	theme(
		panel.grid.major.x = element_line(size = 0.5),
		panel.grid.minor.x = element_line(size = 0.2, linetype = "longdash"),
		axis.text = element_text(size = 12),
		axis.title = element_text(size = 14),
		axis.ticks.y = element_blank(),
		axis.ticks.length = unit(0.6, "lines"),
		plot.title = element_text(hjust = 0.5, size = 16))



# ------------------------------------------------ #
# Identified 5 Ori-Ter symmetrical zones. Boundaries assigned by eye
zoneBounds	<- data.frame(boundary = c(0.045, 0.130, 0.230, 0.375, 0.625, 0.770, 0.870, 0.955, 1.045, 1.130, 1.230))
boundCol	<- wes_palette("Darjeeling1")[1]

# Plot the prelim boundaries
prelimBoundary_plot	<- linearXSpecies_verHGT_dens_plot +
	geom_vline(xintercept = zoneBounds$boundary, color = boundCol) +
	theme(
		panel.grid.major.x = element_line(size = 0.5),
		panel.grid.minor.x = element_blank())

# ------------------------------------------------ #
# The padding up and downstream of a boundary - the transition zone. The total zone is then 45kb @ padding of 0.0075 for a 3mb genome
padding_size	<- 0.0075

# Produce of a set of transition zone coordinates
boundaryPaddding_list	<- lapply(zoneBounds$boundary, function(position) {
	# Padding around transition boundary
	pad_up		<- position - padding_size
	pad_down	<- position + padding_size

	# Return df
	out_df		<- data.frame(boundary = position, boundaryMin = pad_up, boundaryMax = pad_down, stringsAsFactors = FALSE)
	return(out_df)
})
zoneBounds_padded	<- bind_rows(boundaryPaddding_list)


# Plot the prelim boundaries with the transition zone
prelimBoundary_padded_plot	<- linearXSpecies_verHGT_dens_plot +
	# Boundary lines
	geom_vline(xintercept = zoneBounds_padded$boundary, color = boundCol, linetype = "dashed") +
	# Boundary padding
	geom_rect(data = zoneBounds_padded, aes(xmin = boundaryMin, xmax = boundaryMax, ymin = -Inf, ymax = Inf), fill = alpha(boundCol, 0.4), inherit.aes = FALSE) +
	theme(
		panel.grid.major.x = element_line(size = 0.5),
		panel.grid.minor.x = element_blank())


# ------------------------------------------------ #
# Produce of a set of enrichment zones (between boundaries)
zoneRange_list	<- lapply(zoneBounds$boundary, function(position) {

	# Find the next smallest position
	lower_boundaries	<- zoneBounds$boundary[which(zoneBounds$boundary < position)]
	if (length(lower_boundaries) != 0) lower_bound <- max(lower_boundaries) else lower_bound <- 0

	# Zone start and end (padding ignored)
	zoneMin <- lower_bound
	zoneMax	<- position
	
	# Zone start and end (padding included)
	if (identical(lower_bound, 0)) zoneMin_pad	<- zoneMin else zoneMin_pad <- zoneMin + padding_size
	zoneMax_pad	<- zoneMax - padding_size

	# Return df
	out_df		<- data.frame(boundary = position, zoneMin = zoneMin, zoneMax = zoneMax, zoneMin_pad = zoneMin_pad, zoneMax_pad = zoneMax_pad, stringsAsFactors = FALSE)
	return(out_df)
})
zoneRange_df	<- bind_rows(zoneRange_list)
zoneRange_df$zoneName	<- c("oriVer", "oriHGT", "midOther", "flankVer", "terHGT", "flankVer", "midOther", "oriHGT", "oriVer", "oriHGT", "midOther")
zoneRange_df$zoneType	<- c("Ver", "HGT", "Other", "Ver", "HGT", "Ver", "Other", "HGT", "Ver", "HGT", "Other")

zoneCol_df		<- data.frame(
	zoneType = c("HGT", "Ver", "Other"),
	zoneCol = alpha(c(wes_palette("Darjeeling1")[2:3], wes_palette("IsleofDogs1")[6]), 0.2),
	stringsAsFactors = FALSE)
zoneRange_df	<- left_join(zoneRange_df, zoneCol_df, by = "zoneType")

# Finall join the boundary df with the zone df
zoneBoundRange_df	<- left_join(zoneBounds_padded, zoneRange_df, by = "boundary")

# Plot the prelim boundaries with the enrichment zones
prelimBoundary_withZones_plot	<- linearXSpecies_verHGT_dens_plot +
	scale_x_continuous(
		expand = expand_scale(mult = c(0.025, 0.025)),
		name = "Normalized genome position",
		limits = c(0, 1.25),
		breaks = seq(0, 1, by = 0.5),
		labels = c("Origin", "Terminus", "Origin"),
		minor_breaks = seq(0, 1.25, by = ((1 / 360) * 30)),
		sec.axis = dup_axis(
			name = NULL,
			breaks = rowMeans(zoneBoundRange_df[c("zoneMin", "zoneMax")]),
			labels = zoneBoundRange_df$zoneName)
	) +
	# Boundary lines
	geom_vline(xintercept = zoneBoundRange_df$boundary, color = boundCol, linetype = "dashed") +
	# Boundary padding
	geom_rect(data = zoneBoundRange_df, aes(xmin = boundaryMin, xmax = boundaryMax, ymin = -Inf, ymax = Inf), fill = alpha(boundCol, 0.4), inherit.aes = FALSE) +
	# Zone colouring
	geom_rect(data = zoneBoundRange_df, aes(xmin = zoneMin_pad, xmax = zoneMax_pad, ymin = -Inf, ymax = Inf), fill = zoneBoundRange_df$zoneCol, inherit.aes = FALSE) +
	theme(
		panel.grid.major.x = element_line(size = 0.5),
		panel.grid.minor.x = element_blank())


# ------------------------------------------------ #

# Create a truncated zone df - where the final zone ends at 1 (represent whole genome)
finalBoundary		<- subset(zoneBoundRange_df, zoneMin < 1 & zoneMax > 1, select = boundary, drop = TRUE)
zoneRangeToOri_df	<- subset(zoneBoundRange_df, boundary <= finalBoundary)
zoneRangeToOri_df[nrow(zoneRangeToOri_df), c("boundary", "boundaryMin", "boundaryMax", "zoneMax", "zoneMax_pad")] <- c(rep(NA, 3), rep(1, 2))

# Create a truncated zone df - where the final zone ends at 0.5 (for distToOri representations)
finalBoundary		<- subset(zoneBoundRange_df, zoneMin < 0.5 & zoneMax > 0.5, select = boundary, drop = TRUE)
zoneRangeToTer_df	<- subset(zoneBoundRange_df, boundary <= finalBoundary)
zoneRangeToTer_df[nrow(zoneRangeToTer_df), c("boundary", "boundaryMin", "boundaryMax", "zoneMax", "zoneMax_pad")] <- c(rep(NA, 3), rep(0.5, 2))


toOriDensity_withZones_plot	<- ggplot(data = normalPos_df, aes(x = distToOri, color = type)) +
	scale_x_continuous(
		expand = expand_scale(mult = c(0.025, 0.025)),
		name = "Normalized genome position",
		limits = c(0, 0.5),
		breaks = seq(0, 0.5, by = 0.5),
		labels = c("Origin", "Terminus"),
		minor_breaks = seq(0, 0.5, by = ((1 / 360) * 30)),
		sec.axis = dup_axis(
			name = NULL,
			breaks = rowMeans(zoneRangeToTer_df[c("zoneMin", "zoneMax")]),
			labels = zoneRangeToTer_df$zoneName)
	) +	
	scale_y_continuous(name = "Gene Enrichment") +
	# Plot density lines
	stat_density(geom = "line", position = "identity", n = 2^12, adjust = 1/10, size = 1) +
	# Boundary lines
	geom_vline(xintercept = zoneRangeToTer_df$boundary, color = boundCol, linetype = "dashed") +
	# Boundary padding
	geom_rect(data = zoneRangeToTer_df, aes(xmin = boundaryMin, xmax = boundaryMax, ymin = -Inf, ymax = Inf), fill = alpha(boundCol, 0.4), inherit.aes = FALSE) +
	# Zone colouring
	geom_rect(data = zoneRangeToTer_df, aes(xmin = zoneMin_pad, xmax = zoneMax_pad, ymin = -Inf, ymax = Inf), fill = zoneRangeToTer_df$zoneCol, inherit.aes = FALSE) +
	# Coloring, title, themes etc..
	scale_color_manual(values = all_HGT_Ver_cols) +
	ggtitle("Exploratory zone boundary analysis") +
	theme_classic() +
	theme(
		panel.grid.major.x = element_line(size = 0.5),
		panel.grid.minor.x = element_blank(),
		axis.text = element_text(size = 12),
		axis.title = element_text(size = 14),
		axis.ticks.y = element_blank(),
		axis.ticks.length = unit(0.6, "lines"),
		plot.title = element_text(hjust = 0.5, size = 16)
	)


# ------------------------------------------------------------------------------------- #



AG_noPlasmid_path	<- file.path(genome_path, "AG_genome_noPlasmid_gbffs")
AG_noPlasmid_gbk_df	<- data.frame(fileNames = dir(path = AG_noPlasmid_path, pattern = "*.gbk$"), stringsAsFactors = FALSE)
AG_noPlasmid_gbk_df$acc_ass	<- str_replace(AG_noPlasmid_gbk_df$fileNames, pattern="_genomic.gbk", replacement = "")

# Open the taxid / acc_ass translation table
accAss_taxid_file	<- file.path(genome_path, "Genome_lists", "Acc_ass_taxid_table.tsv")
accAss_taxid_df		<- read.table(file = accAss_taxid_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
names(accAss_taxid_df)	<- c("acc_ass", "taxid")

# Find the taxids for these files
AG_noPlasmid_gbk_df	<- left_join(AG_noPlasmid_gbk_df, accAss_taxid_df, by = "acc_ass")

# Remove the rearrange genomes
AG_noPlasmid_gbk_df	<- subset(AG_noPlasmid_gbk_df, !taxid %in% outlierTaxid)

# Create an output folder to hold the tRNA fasta
tRNA_output_path	<- file.path(giProcess_path, "tRNA_fasta")
if (!dir.exists(tRNA_output_path)) dir.create(tRNA_output_path)

# Extract the tRNA positions from the genbank file
invisible(lapply(1:nrow(AG_noPlasmid_gbk_df), function(index) {

	# Entry and output file name
	AG_genbankEntry	<- AG_noPlasmid_gbk_df[index,]
	outputName		<- paste0(AG_genbankEntry$acc_ass, "_tRNA.fasta")

	# See the python script for usage instructions
	system2("genbank_to_fasta.py",
		args = paste0(
			"-i ", file.path(AG_noPlasmid_path, AG_genbankEntry$fileNames),
			" -o ", file.path(tRNA_output_path, outputName),
			" -m 'genbank' -s 'nt' -f 'tRNA' -d 'spacepipe' -q 'protein_id,locus_tag,gene,product,location' &>/dev/null"
		)
	)
}))

















# ------------------------------------------------ #

# How much of the genome are the boundary zones?
transZoneTotal_length	<- sum(zoneBoundary_df$boundaryMax - zoneBoundary_df$boundaryMin, na.rm = TRUE)

HTgenesInTransZones_list	<- lapply(1:nrow(zoneBoundary_df), function(boundaryIndex) {

	boundaryRow	<- zoneBoundary_df[boundaryIndex,]
	boundaryVal	<- boundaryRow$boundary

	transZoneMin	<- boundaryRow$boundaryMin
	transZoneMax	<- boundaryRow$boundaryMax

	if (is.na(boundaryVal)) return(NA)

	HGTgenesInTransZone <- subset(perTypeData$lHGT$'4'$allPosData, relGeneStart >= transZoneMin & relGeneEnd >= transZoneMin & relGeneStart <= transZoneMax & relGeneEnd <= transZoneMax & !taxid %in% outlierTaxid, select = -c(CircStart, CircEnd))

	HGTgenesInTransZone$transZoneID	<- boundaryRow$boundIndex

	return(HGTgenesInTransZone)
})

HTgenesInTransZones_list	<- HTgenesInTransZones_list[!is.na(HTgenesInTransZones_list)]
HTgenesInTransZones_df	<- bind_rows(HTgenesInTransZones_list)



HTgenesInZones_list	<- lapply(1:nrow(zoneBoundary_df), function(boundaryIndex) {

	boundaryRow	<- zoneBoundary_df[boundaryIndex,]

	zoneMin	<- boundaryRow$zoneMin
	zoneMax	<- boundaryRow$zoneMax

	HGTgenesInZone <- subset(perTypeData$lHGT$'4'$allPosData, relGeneStart >= zoneMin & relGeneEnd >= zoneMin & relGeneStart <= zoneMax & relGeneEnd <= zoneMax & !taxid %in% outlierTaxid, select = -c(CircStart, CircEnd))

	HGTgenesInZone$zoneID	<- boundaryRow$zoneType

	return(HGTgenesInZone)
})

HTgenesInZones_list	<- HTgenesInZones_list[!is.na(HTgenesInZones_list)]
HTgenesInZones_df	<- bind_rows(HTgenesInZones_list)







## Non-rearranged HT genes
HTgenes_noOutlier	<- subset(perTypeData$lHGT$`4`$allPosData, !taxid %in% outlierTaxid)

# Total HT genes in 23 non-rearranged genomes == 9034. 1838 orthGroups. 1238 recent HTgenes (13.4%)
nrow(HTgenes_noOutlier)
length(unique(HTgenes_noOutlier$orthGroup))
length(which(HTgenes_noOutlier$Subgroup == TRUE))

# Total HT genes within transition zones in 23 non-rearranged genomes == 983 (10.8% of all HT genes). 415 orthGroups (22.5% of all orthGroups). 115 recent HTgenes (27.7%)
nrow(HTgenesInTransZones_df)
length(unique(HTgenesInTransZones_df$orthGroup))
length(which(HTgenesInTransZones_df$Subgroup == TRUE))


ggplot(data = subset(HTgenes_noOutlier, Subgroup == TRUE), aes(x = relGeneStart)) +
	scale_x_continuous(
		name = "Normalized genome position",
		limits = c(0, 1),
		breaks = seq(0, 1, by = 0.5),
		labels = c("Origin", "Terminus", "Origin"),
		minor_breaks = seq(0, 1.25, by = 0.05)) +
	geom_density(adjust = 1/20) +
	# Boundary lines
	geom_vline(xintercept = zoneBoundary_df$boundary, color = boundaryCol, linetype = "dashed") +
	# Boundary padding
	geom_rect(data = zoneBoundary_df, aes(xmin = boundaryMin, xmax = boundaryMax, ymin = -Inf, ymax = Inf), fill = alpha(boundaryCol, 0.4), inherit.aes = FALSE) +
	# Zone colouring
	geom_rect(data = zoneBoundary_df, aes(xmin = zoneMin, xmax = zoneMax, ymin = -Inf, ymax = Inf), fill = zone_cols[1:9], inherit.aes = FALSE) +
	theme_classic() +
	theme(
		panel.grid.major.x = element_line(size = 0.5),
		panel.grid.minor.x = element_line(size = 0.0, linetype = "longdash"),
		axis.text = element_text(size = 12),
		axis.title = element_text(size = 14))



## Quck COG
transZoneCOG_df	<- as.data.frame(sort(table(unlist(HTgenesInTransZones_df$COGcat)), decreasing = TRUE))
names(transZoneCOG_df)	<- c("COG", "transZone")

allHGTCOGs_df	<- as.data.frame(sort(table(unlist(HTgenes_noOutlier$COGcat)), decreasing = TRUE))
names(allHGTCOGs_df)	<- c("COG", "all_HTgenes")


ggplot(data = subset(HTgenes_noOutlier, COGcat == "M"), aes(x = relGeneStart)) +
	scale_x_continuous(
		name = "Normalized genome position",
		limits = c(0, 1),
		breaks = seq(0, 1, by = 0.5),
		labels = c("Origin", "Terminus", "Origin"),
		minor_breaks = seq(0, 1.25, by = 0.05)) +
	geom_density(n = 2^12, adjust = 1/50) +
	# Boundary lines
	geom_vline(xintercept = zoneBoundary_df$boundary, color = boundaryCol, linetype = "dashed") +
	# Boundary padding
	geom_rect(data = zoneBoundary_df, aes(xmin = boundaryMin, xmax = boundaryMax, ymin = -Inf, ymax = Inf), fill = alpha(boundaryCol, 0.4), inherit.aes = FALSE) +
	# Zone colouring
	geom_rect(data = zoneBoundary_df, aes(xmin = zoneMin, xmax = zoneMax, ymin = -Inf, ymax = Inf), fill = zone_cols[1:9], inherit.aes = FALSE) +
	theme_classic() +
	theme(
		panel.grid.major.x = element_line(size = 0.5),
		panel.grid.minor.x = element_line(size = 0.0, linetype = "longdash"),
		axis.text = element_text(size = 12),
		axis.title = element_text(size = 14))





nearOriHGT_M_genes		<- subset(HTgenesInZones_df, COGcat == "M")
nearOriHGT_M_protIDs	<- nearOriHGT_M_genes$protID

nearOriHGT_M_prod_tbl	<- dbSendQuery(dbConn, 'SELECT * FROM t1 WHERE protID = :protIDs')
dbBind(nearOriHGT_M_prod_tbl, param = list(protIDs = nearOriHGT_M_protIDs))
nearOriHGT_M_prod_df	<- dbFetch(nearOriHGT_M_prod_tbl)
dbClearResult(nearOriHGT_M_prod_tbl)

nearOriHGT_M_groups		<- unique(nearOriHGT_M_prod_df$OrthGroup)
# Are all orthologs of these genes in the same region in Geobacillus spp?
nearOriHGT_M_allGenes	<- subset(perTypeData$lHGT$'4'$allPosData, orthGroup %in% nearOriHGT_M_groups & !taxid %in% outlierTaxid, select = -c(CircStart, CircEnd))

ggplot(data = nearOriHGT_M_allGenes, aes(x = relGeneStart)) +
	scale_x_continuous(
		name = "Normalized genome position",
		limits = c(0, 1),
		breaks = seq(0, 1, by = 0.5),
		labels = c("Origin", "Terminus", "Origin"),
		minor_breaks = seq(0, 1.25, by = 0.05)) +
	geom_density(n = 2^12, adjust = 1/50) +
	# Boundary lines
	geom_vline(xintercept = zoneBoundary_df$boundary, color = boundaryCol, linetype = "dashed") +
	# Boundary padding
	geom_rect(data = zoneBoundary_df, aes(xmin = boundaryMin, xmax = boundaryMax, ymin = -Inf, ymax = Inf), fill = alpha(boundaryCol, 0.4), inherit.aes = FALSE) +
	# Zone colouring
	geom_rect(data = zoneBoundary_df, aes(xmin = zoneMin, xmax = zoneMax, ymin = -Inf, ymax = Inf), fill = zone_cols[1:9], inherit.aes = FALSE) +
	theme_classic() +
	theme(
		panel.grid.major.x = element_line(size = 0.5),
		panel.grid.minor.x = element_line(size = 0.0, linetype = "longdash"),
		axis.text = element_text(size = 12),
		axis.title = element_text(size = 14))




