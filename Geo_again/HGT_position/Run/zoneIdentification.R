#!/usr/bin/env Rscript

# Load master variables and HGT position functions
invisible(sapply(HGTPos.all, source, .GlobalEnv))

# Any libraries explicitly used in the script
require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load("RSQLite", "dplyr","ggplot2", "ggpubr", "wesanderson", "gridExtra")


# ------------------------------------------------------------------------------------- #
# Open All_prot database
dbConn			<- dbConnect(RSQLite::SQLite(), allProtDB_path)

# ------------------------------------------------ #

dataTypes_trunc		<- c("All", "lHGT", "Ver")

extendedPos_data <- lapply(dataTypes_trunc, function(dataType) {

	byAge	<- FALSE
	if (identical(dataType, "All")) {
		penalty		<- NULL
	} else if (identical(dataType, "Ver")) {
		penalty		<- verPenalty
	} else {
		penalty		<- hgtPenalty
		if (identical(dataType, "Old") || identical(dataType, "Recent")) {
			byAge		<- TRUE
			subgroup	<- ifelse(identical(dataType, "Old"), FALSE, TRUE)
		}
	}

	# Get the data per Type
	if (is.null(penalty)) {
		data_trunc	<- subset(perTypeData$All$allPosData, select = -c(CircStart, CircEnd))
	} else if (!byAge) {
		data_trunc	<- subset(perTypeData[[dataType]][[penalty]]$allPosData, select = -c(CircStart, CircEnd))
	}

	# Remove the two outlier genomes
	data_trunc		<- subset(data_trunc, !taxid %in% outlierTaxid)

	data_trunc$type	<- dataType

	firstQ		<- subset(data_trunc, relGeneStart >= 0 & relGeneStart <= 0.25)
	firstQ$relGeneStart	<- firstQ$relGeneStart + 1
	firstQ$relGeneEnd	<- firstQ$relGeneEnd + 1
	extended	<- bind_rows(list(data_trunc, firstQ))

	return(extended)
})

# Convert to df and factor the dataType for colouring
extendedPos_df		<- bind_rows(extendedPos_data)
extendedPos_df$type	<- factor(extendedPos_df$type, levels = dataTypes_trunc)


# ------------------------------------------------ #


# Transition zone boundaries - assigned by eye
boundary_vector	<- c(0.045, 0.955, 1.045, 0.130, 0.870, 1.130, 0.230, 0.770, 1.230, 0.375, 0.625)
# The padding up and downstream of a boundary - the transition zone. The total zone is then 45kb @ padding of 0.0075 for a 3mb genome
padding_size	<- 0.0075

# Produce of a set of transition zone coordinates
boundaryPad_list	<- lapply(boundary_vector, function(position) {
	# Padding around transition boundary
	pad_up		<- position - padding_size
	pad_down	<- position + padding_size

	# Return df
	out_df		<- data.frame(boundary = position, boundaryMin = pad_up, boundaryMax = pad_down, ymin = -Inf, ymax = Inf, stringsAsFactors = FALSE)
	return(out_df)
})
boundaryPad_df	<- bind_rows(boundaryPad_list)

# Produce of a set of enrichment zones (between transition zones)
zoneRange_list	<- lapply(sort(boundary_vector), function(position) {
	# Find the next smallest position
	lower_boundaries	<- boundary_vector[which(boundary_vector < position)]
	if (length(lower_boundaries) != 0) lower_bound <- max(lower_boundaries) else lower_bound <- 0
	
	# If it's the first zone from the origin, the lower bound is 0 (unpadded)
	if (identical(lower_bound, 0)) zoneMin	<- lower_bound else zoneMin <- lower_bound + padding_size
	zoneMax	<- position - padding_size

	# Return df
	out_df		<- data.frame(boundary = position, zoneMin = zoneMin, zoneMax = zoneMax, ymin = -Inf, ymax = Inf, stringsAsFactors = FALSE)
	return(out_df)
})
zoneRange_df	<- bind_rows(zoneRange_list)

# Combine into a single dataframe
zoneBoundaryExt_df	<- cbind(merge(boundaryPad_df[-4:-5], zoneRange_df[-4:-5], by = "boundary"), zoneRange_df[,4:5])

## Create a truncated zone df - where the final zone ends at 1 (represent whole genome)
finalBoundary		<- subset(zoneBoundaryExt_df, zoneMin < 1 & zoneMax > 1, select = boundary, drop = TRUE)
# Drop the repeat zones and ymin/ymax needed for plotting
zoneBoundary_df		<- subset(zoneBoundaryExt_df, boundary <= finalBoundary, select = -c(ymin, ymax))
# Adjust final zone to end at 1
zoneBoundary_df[which(zoneBoundary_df$boundary == finalBoundary),] <- c(NA, NA, NA, zoneBoundary_df$zoneMin[which(zoneBoundary_df$boundary == finalBoundary)], 1)
# Add the type of zone
zoneBoundary_df$zoneType	<- c("oriVer", "oriHGT", "midOther", "flankVer", "terHGT", "flankVer", "midOther", "oriHGT", "oriVer")
zoneBoundary_df$boundIndex	<- c("a1", "b1", "c1", "d1", "d2", "c2", "b2", "a2", NA)

# ------------------------------------------------ #
# Plot the density curves of HGT and vertical genes and boundaries

# Density curve colours
density_cols	<- c(alpha(wes_palette("Chevalier1")[4], 1), wes_palette("Darjeeling1")[2:3])

# Set up zone colours
hgtZoneCol		<- alpha(wes_palette("Darjeeling1")[2], 0.2)
verZoneCol		<- alpha(wes_palette("Darjeeling1")[3], 0.2)
othZoneCol		<- alpha(wes_palette("IsleofDogs1")[6], 0.2)
zone_cols		<- c(verZoneCol, hgtZoneCol, othZoneCol, verZoneCol, hgtZoneCol, verZoneCol, othZoneCol, hgtZoneCol, verZoneCol, hgtZoneCol, othZoneCol)

# Boundary color
boundaryCol		<- wes_palette("Darjeeling1")[1]

# Plot
splitRelGenome_zones <- ggplot(data = plot_df, aes(x = relGeneStart, color = type)) +
	scale_x_continuous(
		name = "Normalized genome position",
		limits = c(0, 1.25),
		breaks = seq(0, 1, by = 0.5),
		labels = c("Origin", "Terminus", "Origin"),
		minor_breaks = seq(0, 1.25, by = 0.05)) +
	scale_y_continuous(name = "Gene Density") +
	stat_density(geom = "line", position = "identity", n = 2^12, adjust = 1/12, size = 1) +
	scale_color_manual(values = plot_colors) +
	# Boundary lines
	geom_vline(xintercept = zoneBoundaryExt_df$boundary, color = boundaryCol, linetype = "dashed") +
	# Boundary padding
	geom_rect(data = zoneBoundaryExt_df, aes(xmin = boundaryMin, xmax = boundaryMax, ymin = ymin, ymax = ymax), fill = alpha(boundaryCol, 0.4), inherit.aes = FALSE) +
	# Zone colouring
	geom_rect(data = zoneBoundaryExt_df, aes(xmin = zoneMin, xmax = zoneMax, ymin = ymin, ymax = ymax), fill = zone_cols, inherit.aes = FALSE) +
	theme_classic() +
	theme(
		panel.grid.major.x = element_line(size = 0.5),
		panel.grid.minor.x = element_line(size = 0.0, linetype = "longdash"),
		axis.text = element_text(size = 12),
		axis.title = element_text(size = 14))

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




