#!/usr/bin/env Rscript

# Load master variables and HGT position functions
invisible(sapply(HGTPos.all, source, .GlobalEnv))

# Any libraries explicitly used in the script
require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load("RSQLite", "dplyr","ggplot2", "ggtree", "wesanderson", "Biostrings", "ggpubr", "Hmisc", "data.table", "broom")


# ------------------------------------------------------------------------------------- #
# Read in data

message("\nReading in data...", appendLF = FALSE)

# General position data
perTypeData			<- readRDS(file.path(positionData_path, "AG_perTypeData.rds"))

# Per-type COG data
perTypeCOG_data		<- readRDS(file.path(positionData_path, "AG_perTypeCOGData.rds"))

# Subgroup data
subgroupData		<- readRDS(file.path(positionData_path, "AG_subgroupData.rds"))

# Genomic Island data
byType_withGI_data	<- readRDS(file.path(positionData_path, "AG_perTypeGIData.rds"))

# Zone list
zoneBoundaryList	<- readRDS(file.path(positionData_path, "AG_zoneBoundaries.rds"))

message("\rReading in data... done\n")

# ------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------- #

# Important variables
contPalette		<- colorRampPalette(wes_palette("Zissou1"))
binomial_list	<- unique(perTypeData$All$allPosData$binomial)

binNumber		<- 200
minGeneNum		<- 3

# Open All_prot database
conn			<- dbConnect(RSQLite::SQLite(), allProtDB_path)

# ------------------------------------------------------------------------------------- #
# Plotting

# Output path for figures
GC_TreesFig_path	<- file.path(figureOutput_path, "GC_content", "Trees")
GC_SpeciesFig_path	<- file.path(figureOutput_path, "GC_content", "bySpecies")
GC_byBinFig_path	<- file.path(figureOutput_path, "GC_content", "byBin")

if (!dir.exists(GC_TreesFig_path)) dir.create(GC_TreesFig_path, recursive = TRUE)
if (!dir.exists(GC_SpeciesFig_path)) dir.create(GC_SpeciesFig_path)
if (!dir.exists(GC_byBinFig_path)) dir.create(GC_byBinFig_path)

# Quartz plotting options common to this script
quartz.options(canvas = "white", bg = "white")


# ------------------------------------------------------------------------------------- #
# ---	Prepare phylogenies for plotting output - need tip order for plotting later	--- #
# ------------------------------------------------------------------------------------- #

# Get the binomial tree from the subgroup data set
AG_binomTime_tree	<- subgroupData$AG_binomTime_tree

# Plot with and without branch.lenght to get phylo and clado trees
binomTree_phylo <- ggtree(tr = AG_binomTime_tree, ladderize = FALSE, col =  axisCol, size = 1.5) +
	geom_tiplab(angle = -90, hjust = -0.05, col = textCol) +
	scale_x_reverse(limits = c(0.5, 0)) +
	coord_flip() +
	lightTheme +
	theme(
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank()
	)

binomTree_clado <- ggtree(tr = AG_binomTime_tree, branch.length = "none", ladderize = FALSE, col =  axisCol, size = 1.5) +
	geom_tiplab(angle = -90, hjust = -0.05, col = textCol) +
	scale_x_reverse(limits = c(20, 0)) +
	coord_flip() +
	lightTheme +
	theme(
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank()
	)

# Prepare tree as dataframe to extract tip order
treeAsDf	<- subset(fortify(AG_binomTime_tree, ladderize = FALSE), isTip)
tipOrder	<- with(treeAsDf, label[order(y, decreasing = T)])

# Write trees out to PDF for combination plots
quartz(width = 22, height = 12)
print(binomTree_phylo)
quartz.save(file = file.path(GC_TreesFig_path, "BinomialTree_phylo.pdf"), type = "pdf", dpi = 300)
invisible(dev.off())

quartz(width = 16, height = 8)
print(binomTree_clado)
quartz.save(file = file.path(GC_TreesFig_path, "BinomialTree_clado.pdf"), type = "pdf", dpi = 300)
invisible(dev.off())


# ------------------------------------------------------------------------------------- #
# ---	Prepare GC-content data structure - this is done over N (100) genomic bins	--- #
# ------------------------------------------------------------------------------------- #

message("Calculating GC content profiles for all dataTypes in parallel")
byTypebySpeciesGC_data	<- mclapply(dataTypes_withAge, function(dataType) {

	# Set up variables depending on the dataType
	byAge		<- FALSE
	subgroup	<- NA 
	geneCutOff	<- 1

	if (identical(dataType, "All")) {
		penalty		<- NA
	} else if (identical(dataType, "Ver")) {
		penalty		<- "3"
		geneCutOff	<- minGeneNum
	} else {
		penalty		<- "4"
		if (identical(dataType, "Old") || identical(dataType, "Recent")) {
			byAge		<- TRUE
			subgroup	<- ifelse(identical(dataType, "Old"), FALSE, TRUE)
		}
	}

	# Iterate over species, each genome divided into N equally sized bins
	bySpeciesGC_list	<- lapply(binomial_list, function(species) {

		# For each bin, isolate genes in that bin and calculate their GC content
		byBinGC_list	<- lapply(seq(1:binNumber), function(genome_bin) {

			# Define bin position
			bin_start		<- (genome_bin - 1) / binNumber
			bin_end			<- genome_bin / binNumber

			# Isolate the data we need for this dataType
			if (is.na(penalty)) {
				subsetData		<- perTypeData$All$allPosData
			} else if (!byAge) {
				subsetData		<- perTypeData[[dataType]][[penalty]]$allPosData
			} else {
				subsetData		<- subset(perTypeData$lHGT[[penalty]]$allPosData, Subgroup == subgroup)
			}

			# Further subset by position
			subsetByPos			<- subset(subsetData, binomial == species & relGeneStart >= bin_start & relGeneStart <= bin_end)
			subsetByPos			<- subsetByPos[order(subsetByPos$relGeneStart),]

			# If the number of genes of this dataType in this Bin is 0 - return
			if (length(subsetByPos$protID) < geneCutOff) {
				return(data.frame(BinIndex = bin_end, GC_content = NA, stringsAsFactors = FALSE))
			}

			# Get the corresponding nucleotide sequences from the Sqlite database
			nuclSeq_tbl		<- dbSendQuery(conn, 'SELECT NuclSeq FROM t1 WHERE protID = :protIDs AND locus = :locTags AND binomial = :species')
			dbBind(nuclSeq_tbl, param = list(protIDs = subsetByPos$protID, locTags = subsetByPos$locusTag, species = rep(species, nrow(subsetByPos))))
			nuclSeq_df		<- dbFetch(nuclSeq_tbl)
			dbClearResult(nuclSeq_tbl)

			# Convert to DNAStringSet
			nuclSeq_stringSet	<- DNAStringSet(x = nuclSeq_df$NuclSeq)

			# For lHGT genes - calculate overall GC content
			GeneGC_cont		<- rowSums(subset(alphabetFrequency(nuclSeq_stringSet, as.prob = TRUE), select = 2:3, drop = FALSE))
			
			out_df			<- data.frame(BinIndex = bin_end, GC_content = GeneGC_cont, locusTag = subsetByPos$locusTag, stringsAsFactors = FALSE)

			return(out_df)
		})

		# Combine the per-Bin data
		byBinGC_df			<- bind_rows(byBinGC_list)
		byBinGC_df$Species	<- species

		# Check for the "All" gene groups that every gin makes it into the bin (and not more than once into any bins)
		# if (nrow(byBinGC_df[!is.na(byBinGC_df$binLocTags),]) != nrow(subset(perTypeData$All$allPosData, binomial == species))) message(species)

		return(byBinGC_df)
	})
	# Combine the per-Species data
	bySpeciesGC_df		<- bind_rows(bySpeciesGC_list)
	bySpeciesGC_df$Type	<- dataType

	return(bySpeciesGC_df)
}, mc.cores = 20)

# Rename list after dataTypes
names(byTypebySpeciesGC_data)	<- dataTypes_withAge

# Calculate normalised GC content for cross-species comparison
byTypeSpeciesGCNorm		<- lapply(byTypebySpeciesGC_data, function(GC_byType) {

	bySpeciesNorm_list	<- lapply(binomial_list, function(speciesName) {
		speciesMeanGC	<- mean(subset(byTypebySpeciesGC_data$All, Species == speciesName)$GC_content, na.rm = TRUE)
		bySpecies_data	<- subset(GC_byType, Species == speciesName)
		bySpecies_data$GC_norm	<- bySpecies_data$GC_content - speciesMeanGC
		return(bySpecies_data)
	})
	bySpeciesNorm_df	<- bind_rows(bySpeciesNorm_list)
})

dbDisconnect(conn)



# ------------------------------------------------------------------------------------- #
# ---	Examine overall GC-content trends per-Species - compare Ver/HGT + New/Old	--- #
# ------------------------------------------------------------------------------------- #


# ---	Plot overall (all gene) GC content for each species to compare phylogeny	--- #

GC_all					<- byTypeSpeciesGCNorm$All
GC_all$Species			<- factor(GC_all$Species, levels = rev(tipOrder))

perSpecies_cols	<- colorRampPalette(wes_palette("Darjeeling1"))(length(binomial_list))

GC_all_plot		<- ggplot(data = GC_all, aes(x = Species, y = GC_content, color = Species, fill = Species)) +
	scale_y_continuous(name = "Genic GC content") +
	geom_violin(col = axisCol) +
	scale_fill_manual(values = perSpecies_cols, guide = FALSE) +
	lightTheme +
	theme(
		axis.text.x = element_blank()
	)

quartz(width = 18, height = 8)
print(GC_all_plot)
quartz.save(file = file.path(GC_SpeciesFig_path, "GC_allBySpecies.pdf"), type = "pdf", dpi = 300)
invisible(dev.off())



# ---	Is there a difference in GC content between All lHGTs and Vertical genes?	--- #
GC_VervsHGT			<- bind_rows(list(byTypeSpeciesGCNorm$Ver, byTypeSpeciesGCNorm$lHGT))

# Factor variables and set up colours 
GC_VervsHGT$Type	<- factor(GC_VervsHGT$Type, levels = dataTypes_withAge)
GC_VervsHGT$Species	<- factor(GC_VervsHGT$Species, levels = rev(tipOrder))

# Prepare plot
GC_VervsHGT_plot	<- ggplot(data = GC_VervsHGT, aes(x = Type, y = GC_content, color = Type, fill = Type)) +
	geom_boxplot(color = axisCol) +
	scale_y_continuous(
		limits = c(round_any(min(GC_VervsHGT$GC_content, na.rm = TRUE), 0.1, floor), round_any(max(GC_VervsHGT$GC_content, na.rm = TRUE), 0.1, ceiling)),
		name = "Genic GC content") +
	facet_wrap(~Species, nrow = 1) +
	stat_compare_means(method = "t.test", col = axisCol, label = "p.signif", label.x.npc = "centre", label.y.npc = "bottom", size = 6) +
	scale_fill_manual(values = c(dataTypeCols$HGT, dataTypeCols$Ver)) +
	lightTheme +
	theme(
		panel.spacing = unit(2, "lines"),
		axis.text.x = element_blank(),
		panel.grid.minor = element_blank(),
		panel.grid.major.x = element_blank(),
		axis.ticks.x = element_blank(),
		strip.background = element_blank(),		
		strip.text = element_text(color = textCol, angle = 90, size = 5)		
	)


quartz(width = 24, height = 10)
print(GC_VervsHGT_plot)
quartz.save(file = file.path(GC_SpeciesFig_path, "GC_VervsHGTbySpecies.pdf"), type = "pdf", dpi = 300)
invisible(dev.off())


# ---	Is there a difference in GC content between OLD lHGTs and RECENT lHGTs?		--- #

# NB - we get some "Recent" values for species that don't have Recent ancestor branches
# This is due to secondary AG2AG transfers from other (recent) branches into these taxa

GC_OldvsNew			<- bind_rows(list(byTypeSpeciesGCNorm$Old, byTypeSpeciesGCNorm$Recent))

# Factor variables and set up colours
GC_OldvsNew$Type	<- factor(GC_OldvsNew$Type, levels = dataTypes_withAge)
GC_OldvsNew$Species	<- factor(GC_OldvsNew$Species, levels = rev(tipOrder))

# Here we count the number of observations for old and recent HGTs per species
bySpeciesOvN_count		<- lapply(binomial_list, function(speciesName) {

	byType_count	<- lapply(unique(GC_OldvsNew$Type), function(dataType) {
		byTypeData	<- subset(GC_OldvsNew, Species == speciesName & Type == dataType)
		max			<- max(byTypeData$GC_content, na.rm = TRUE)
		count		<- length(which(!is.na(byTypeData$GC_content)))
		return(data.frame(Type = dataType, Species = speciesName, Count = count, MaxVal = max, stringsAsFactors = FALSE))
	})

	return(bind_rows(byType_count))
})

bySpeciesOvN_count_df	<- bind_rows(bySpeciesOvN_count) 
bySpeciesOvN_count_df$Species	<- factor(bySpeciesOvN_count_df$Species, levels = rev(tipOrder))

# Plot
GC_OldvsNew_plot	<- ggplot(data = GC_OldvsNew, aes(x = Type, y = GC_content, color = Type, fill = Type)) +
	geom_boxplot(color = axisCol) +
	scale_y_continuous(
		limits = c(round_any(min(GC_OldvsNew$GC_content, na.rm = TRUE), 0.1, floor), round_any(max(GC_OldvsNew$GC_content, na.rm = TRUE), 0.1, ceiling)),
		name = "Genic GC content") +
	facet_wrap(~Species, nrow = 1) +
	geom_text(data = bySpeciesOvN_count_df, aes(y = MaxVal, label = Count), vjust = -1) +
	stat_compare_means(method = "t.test", col = axisCol, label = "p.signif", label.x.npc = "centre", label.y.npc = "bottom", size = 6) +
	scale_fill_manual(values = c(dataTypeCols$Old, dataTypeCols$Recent, dataTypeCols$Old)) +
	scale_color_manual(values = c(dataTypeCols$Old, dataTypeCols$Recent)) +
	lightTheme + 
	theme(
		panel.spacing = unit(2, "lines"),
		axis.text.x = element_blank(),
		panel.grid.minor = element_blank(),
		panel.grid.major.x = element_blank(),
		axis.ticks.x = element_blank(),
		strip.background = element_blank(),		
		strip.text = element_text(color = textCol, angle = 90, size = 5)		
	)

# Write plot
quartz(width = 24, height = 10)
print(GC_OldvsNew_plot)
quartz.save(file = file.path(GC_SpeciesFig_path, "GC_OldvsNewbySpecies.pdf"), type = "pdf", dpi = 300)
invisible(dev.off())


# ---	Is there a difference in GC content between OLD lHGTs and Vertical genes?	--- #

GC_VervsOld			<- bind_rows(list(byTypeSpeciesGCNorm$Ver, byTypeSpeciesGCNorm$Old))

# Factor variables and set up colours
GC_VervsOld$Type	<- factor(GC_VervsOld$Type, levels = dataTypes_withAge)
GC_VervsOld$Species	<- factor(GC_VervsOld$Species, levels = rev(tipOrder))

bySpeciesVvO_count		<- lapply(binomial_list, function(speciesName) {

	byType_count	<- lapply(unique(GC_VervsOld$Type), function(dataType) {
		byTypeData	<- subset(GC_VervsOld, Species == speciesName & Type == dataType)
		max			<- max(byTypeData$GC_content, na.rm = TRUE)
		count		<- length(which(!is.na(byTypeData$GC_content)))
		return(data.frame(Type = dataType, Species = speciesName, Count = count, MaxVal = max, stringsAsFactors = FALSE))
	})

	return(bind_rows(byType_count))
})
bySpeciesVvO_count_df	<- bind_rows(bySpeciesVvO_count)
bySpeciesVvO_count_df$Species	<- factor(bySpeciesVvO_count_df$Species, levels = rev(tipOrder))

GC_VervsOld_plot	<- ggplot(data = GC_VervsOld, aes(x = Type, y = GC_content, color = Type, fill = Type)) +
	geom_boxplot(color = axisCol) +
		scale_y_continuous(
		limits = c(round_any(min(GC_VervsOld$GC_content, na.rm = TRUE), 0.1, floor), round_any(max(GC_VervsOld$GC_content, na.rm = TRUE), 0.1, ceiling)),
		name = "Genic GC content") +
	facet_wrap(~Species) +
	geom_text(data = bySpeciesVvO_count_df, aes(y = MaxVal, label = Count), vjust = -1) +
	stat_compare_means(method = "t.test", col = axisCol, label = "p.signif", label.x.npc = "centre", label.y.npc = "bottom", size = 6) +
	scale_fill_manual(values = c(dataTypeCols$Ver, dataTypeCols$Old)) +
	scale_color_manual(values = c(dataTypeCols$Ver, dataTypeCols$Old)) +
	lightTheme + 
	theme(
		panel.spacing = unit(2, "lines"),
		axis.text.x = element_blank(),
		panel.grid.minor = element_blank(),
		panel.grid.major.x = element_blank(),
		axis.ticks.x = element_blank(),
		strip.background = element_blank(),		
		strip.text = element_text(color = textCol, size = 6)		
	)

# Write plot
quartz(width = 16, height = 16)
print(GC_VervsOld_plot)
quartz.save(file = file.path(GC_SpeciesFig_path, "GC_VervsOldbySpecies.pdf"), type = "pdf", dpi = 300)
invisible(dev.off())






# ------------------------------------------------------------------------------------- #
# ---	Examine GC-content across the genome by looking at perBin GC differences	--- #
# ------------------------------------------------------------------------------------- #

# ---	Compare all-Gene GC-content change per-Species across the genome Ori-to-Ori	--- #

# Get the relevant data + factor
GC_allGenes			<- byTypeSpeciesGCNorm$All
GC_allGenes_ext		<- rbind(GC_allGenes, cbind(GC_allGenes[which(GC_allGenes$BinIndex > 0 & GC_allGenes$BinIndex <= 0.25), 1, drop = FALSE] + 1, GC_allGenes[which(GC_allGenes$BinIndex > 0 & GC_allGenes$BinIndex <= 0.25), -1]))
GC_allGenes_ext$Type	<- factor(GC_allGenes_ext$Type, levels = dataTypes_withAge)

# Plot. This plot is not very informative, as there is a lot of intraGenome fluctuations and the two GC-distinct groups complicate viewing
GCbyBin_all_plot	<- ggplot(GC_allGenes_ext, aes(x = BinIndex, y = GC_content, color = Type)) +
	scale_x_continuous(
		expand = expand_scale(mult = c(0.025, 0.025)),
		name = "Normalized genome position",
		limits = c(0, 1.25),
		breaks = seq(0, 1.25, by = 0.5),
		labels = c("Origin", "Terminus", "Origin"),
		minor_breaks = seq(0, 1.25, by = ((1 / 360) * 30)),
		sec.axis = dup_axis(
			name = NULL,
			breaks = rowMeans(zoneBoundaryList$expandedRange[c("zoneMin", "zoneMax")]),
			labels = zoneBoundaryList$expandedRange$zoneName)) +
	scale_y_continuous(name = "Genic GC content") +
	geom_vline(xintercept = zoneBoundaryList$expandedRange$boundary, color = zoneBoundaryList$boundCol, linetype = "dashed") +	
	geom_rect(data = zoneBoundaryList$expandedRange, aes(xmin = zoneMin, xmax = zoneMax, ymin = -Inf, ymax = Inf), fill = zoneBoundaryList$expandedRange$zoneCol_alpha, inherit.aes = FALSE) +
	geom_smooth(
		data = GC_allGenes_ext,
		aes(group = Species),
		method = "loess",
		span = 0.1,
		col = alpha(dataTypeCols$All, 0.8),
		size = 0.2,
		se = FALSE) +
	# Overall smooth line (cross-species)
	scale_color_manual(values = dataTypeCols$All, guide = guide_legend(title = "Gene Type")) +
	scale_fill_manual(values = alpha(dataTypeCols$All, 0.7), guide = FALSE) +
	lightTheme

# Write plot
quartz(width = 22, height = 8)
print(GCbyBin_all_plot)
quartz.save(file = file.path(GC_byBinFig_path, "perSpeciesAll_raw.pdf"), type = "pdf", dpi = 300)
invisible(dev.off())


# ---		Compare normalised GC-content change per-Species across the genome 		--- #

# Plot with normalised per-Bin GC values. Have to use GAM smoothhing method (loess memory fail), k=25 looks like happy medium.
GCbyBin_allNorm_plot	<- ggplot(GC_allGenes_ext, aes(x = BinIndex, y = GC_norm, color = Type)) +
	scale_x_continuous(
		expand = expand_scale(mult = c(0.025, 0.025)),
		name = "Normalized genome position",
		limits = c(0, 1.25),
		breaks = seq(0, 1.25, by = 0.5),
		labels = c("Origin", "Terminus", "Origin"),
		minor_breaks = seq(0, 1.25, by = ((1 / 360) * 30)),
		sec.axis = dup_axis(
			name = NULL,
			breaks = rowMeans(zoneBoundaryList$expandedRange[c("zoneMin", "zoneMax")]),
			labels = zoneBoundaryList$expandedRange$zoneName)) +
	scale_y_continuous(
		name = "GC content") +
	geom_vline(xintercept = zoneBoundaryList$expandedRange$boundary, color = zoneBoundaryList$boundCol, linetype = "dashed") +	
	geom_rect(data = zoneBoundaryList$expandedRange, aes(xmin = zoneMin, xmax = zoneMax, ymin = -Inf, ymax = Inf), fill = zoneBoundaryList$expandedRange$zoneCol_alpha, inherit.aes = FALSE) +
	geom_smooth(
		data = GC_allGenes_ext,
		aes(group = Species),
		method = "loess",
		span = 0.1,
		col = alpha(dataTypeCols$All, 0.4),
		size = 0.2,
		se = FALSE) +
	# Overall smooth line (cross-species)
	geom_smooth(method = "gam", formula = y ~ s(x, k = 25), size = 1.5, se = TRUE, span = 0.1, color = dataTypeCols$All) +
	scale_fill_manual(values = alpha(dataTypeCols$All, 0.7), guide = FALSE) +
	coord_cartesian(ylim = c(-0.1, 0.05)) +
	ggtitle("Normalized GC content of all genes across 23 AG genomes") +
	lightTheme


# Write plot
quartz(width = 22, height = 8)
print(GCbyBin_allNorm_plot)
quartz.save(file = file.path(GC_byBinFig_path, "perSpeciesAll_normalised.pdf"), type = "pdf", dpi = 300)
invisible(dev.off())


# ---		Compare normalised GC-content between Ver and lHGT genes per-Species	--- #
# Get the relevant data
GC_VerVsHGT				<- bind_rows(list(byTypeSpeciesGCNorm$Ver, byTypeSpeciesGCNorm$lHGT))

# Add the first quarter of the genome onto the end to smooth over the origin
GC_VerVsHGT_ext		<- rbind(GC_VerVsHGT, cbind(GC_VerVsHGT[which(GC_VerVsHGT$BinIndex > 0 & GC_VerVsHGT$BinIndex <= 0.25), 1, drop = FALSE] + 1, GC_VerVsHGT[which(GC_VerVsHGT$BinIndex > 0 & GC_VerVsHGT$BinIndex <= 0.25), -1]))


# Factor variables + color
GC_VerVsHGT_ext$Type	<- factor(GC_VerVsHGT_ext$Type, levels = dataTypes_withAge)
GC_VerVsHGT_ext$Species	<- factor(GC_VerVsHGT_ext$Species, levels = rev(tipOrder))

# Plot with normalised per-Bin GC values. Have to use GAM smoothhing method (loess memory fail), k=20 looks most like loess.
GCbyBin_VerHGT_norm_plot	<- ggplot(GC_VerVsHGT_ext, aes(x = BinIndex, y = GC_norm, color = Type), fill = Type) +
	scale_x_continuous(
		expand = expand_scale(mult = c(0.025, 0.025)),
		name = "Normalized genome position",
		limits = c(0, 1.25),
		breaks = seq(0, 1.25, by = 0.5),
		labels = c("Origin", "Terminus", "Origin"),
		minor_breaks = seq(0, 1.25, by = ((1 / 360) * 30)),
		sec.axis = dup_axis(
			name = NULL,
			breaks = rowMeans(zoneBoundaryList$expandedRange[c("zoneMin", "zoneMax")]),
			labels = zoneBoundaryList$expandedRange$zoneName)) +
	scale_y_continuous(name = "Normalised GC content") +
	geom_vline(xintercept = zoneBoundaryList$expandedRange$boundary, color = zoneBoundaryList$boundCol, linetype = "dashed") +	
	geom_rect(data = zoneBoundaryList$expandedRange, aes(xmin = zoneMin, xmax = zoneMax, ymin = -Inf, ymax = Inf), fill = zoneBoundaryList$expandedRange$zoneCol_alpha, inherit.aes = FALSE) +
	geom_smooth(
		data = subset(GC_VerVsHGT_ext, Type == "Ver"),
		aes(group = Species),
		method = "loess",
		span = 0.1,
		col = alpha(dataTypeCols$Ver, 0.6),
		size = 0.2,
		se = FALSE) +
	geom_smooth(
		data = subset(GC_VerVsHGT_ext, Type == "lHGT"),
		aes(group = Species),
		method = "loess",
		span = 0.1,
		col = alpha(dataTypeCols$HGT, 0.6),
		size = 0.2,
		se = FALSE) +
	# Overall smooth line (cross-species)
	# geom_smooth(aes(color = Type, fill = Type), method = "gam", formula = y ~ s(x, k = 20), size = 1.5, se = TRUE, span = 0.05) +
	geom_smooth(aes(color = Type, fill = Type), method = "loess", size = 1.5, se = TRUE, span = 0.1) +
	scale_color_manual(values = c(dataTypeCols$HGT, dataTypeCols$Ver), guide = guide_legend(title = "Gene Type")) +
	scale_fill_manual(values = alpha(c(dataTypeCols$HGT, dataTypeCols$Ver), 0.7), guide = FALSE) +
	lightTheme +
	ggtitle("Normalized GC content of Vertical and HGT genes across 23 AG genomes") +
	theme(
		legend.justification = c(1, 1),
		legend.position = c(0.95, 0.95),
		legend.background = element_rect(fill = "white", color = axisCol)
	)

# Write plot
quartz(width = 22, height = 10)
print(GCbyBin_VerHGT_norm_plot)
quartz.save(file = file.path(GC_byBinFig_path, "perSpeciesVerHGT_normalised.pdf"), type = "pdf", dpi = 300)
invisible(dev.off())




# ---	Compare Vertical v lHGT GC-content per-Species across the genome Ori-to-Ori	--- #
GCbySpecies_VerVsHGT_plot	<- ggplot(GC_VerVsHGT_ext, aes(x = BinIndex, y = GC_content, color = Type)) +
	scale_x_continuous(
		expand = expand_scale(mult = c(0.025, 0.025)),
		name = "Normalized genome position",
		limits = c(0, 1.25),
		breaks = seq(0, 1.25, by = 0.5),
		labels = c("Origin", "Terminus", "Origin"),
		minor_breaks = seq(0, 1.25, by = ((1 / 360) * 30))) +
	scale_y_continuous(name = "GC content") +
	geom_vline(xintercept = zoneBoundaryList$expandedRange$boundary, color = zoneBoundaryList$boundCol, linetype = "dashed") +	
	geom_rect(data = zoneBoundaryList$expandedRange, aes(xmin = zoneMin, xmax = zoneMax, ymin = -Inf, ymax = Inf), fill = rep(zoneBoundaryList$expandedRange$zoneCol_alpha, length(binomial_list)), inherit.aes = FALSE) +
	geom_smooth(
		data = subset(GC_VerVsHGT_ext, Type == "lHGT"),
		aes(group = Species),
		method = "loess",
		span = 0.1,
		color = dataTypeCols$HGT,
		size = 0.7,
		se = FALSE) +
	geom_smooth(
		data = subset(GC_VerVsHGT_ext, Type == "Ver"),
		aes(group = Species),
		method = "loess",
		span = 0.1,
		color = dataTypeCols$Ver,
		size = 0.7,
		se = FALSE) +
	facet_wrap(~Species, scales = "free") +
	lightTheme +
	theme(
		panel.spacing = unit(2, "lines"),
		panel.grid.minor = element_blank(),
		strip.background = element_rect(fill = "transparent", color = axisCol),		
		strip.text = element_text(color = textCol)
	)

# Write plot
quartz(width = 24, height = 14)
print(GCbySpecies_VerVsHGT_plot)
quartz.save(file = file.path(GC_byBinFig_path, "bySpeciesVerHGT_facet.pdf"), type = "pdf", dpi = 300)
invisible(dev.off())




# ---	Correlate Vertical & lHGT GC-content by Species over the genome Ori-to-Ori	--- #
# Get the relevant data - don't use the one above (normalised to 100% x-axis)
GC_VerVsHGT				<- bind_rows(list(byTypeSpeciesGCNorm$Ver, byTypeSpeciesGCNorm$lHGT))

# Choose correlation method
corrMethod	<- "pearson"

meanBySpecies <- GC_VerVsHGT %>%
	group_by(BinIndex, Species, Type) %>%
	dplyr::summarise(meanGC = mean(GC_content, na.rm = TRUE)) %>%
	spread(Type, meanGC)

corBySpecies <- meanBySpecies %>%
	group_by(Species) %>%
	do(tidy(rcorr(x = .$lHGT, y = .$Ver, type = corrMethod))) %>%
	mutate(pLabel = ifelse(p.value <= 0.05, paste0(" ~italic(p) == ", signif(p.value, digits = 2)), paste0(" ~NS"))) %>%
	mutate(corType = ifelse(identical(corrMethod, "pearson"), paste0("italic(r) == "), paste0("italic(rs) == "))) %>%
	mutate(fullLabel = paste0(corType, signif(estimate, 2), pLabel))


VervHGT_GCCorr_plot	<- ggplot(data = meanBySpecies, aes(x = Ver, y = lHGT, group = Species)) + 
	scale_y_continuous(name = "lHGT GC Content per-Bin") +
	scale_x_continuous(name = "Vertical GC Content per-Bin") +
	geom_point(size = 1, color = axisCol) +
	geom_smooth(method = "lm", color = "red") +
	facet_wrap(~Species) +
	# With scales free, use the pos.x and pos.y values
	geom_text(data = corBySpecies, aes(x = 0.45, y = 0.55, label = fullLabel), parse = TRUE, size = 4, color = "red") +
	lightTheme +
	theme(
		panel.spacing = unit(2, "lines"),
		panel.grid.minor = element_blank(),
		strip.background = element_rect(fill = "transparent", color = axisCol),		
		strip.text = element_text(color = textCol)
	)

# Write plot
quartz(width = 24, height = 14)
print(VervHGT_GCCorr_plot)
quartz.save(file = file.path(GC_byBinFig_path, "correlateBySpecies_VerHGT_facet.pdf"), type = "pdf", dpi = 300)
invisible(dev.off())

# ------------------------------------------------------------------------------------- #
# Does the genomic context of HGT gain location affect either the vertical genes or HGT genes gained?
# Compare the GC content of recent HGT gains across genomes (normalised GC), in bins where there are vertical
# genes and bins where there are no vertical genes

# Need Vertical and recent HGT data
GC_VerVsNew		<- bind_rows(list(byTypeSpeciesGCNorm$Ver, byTypeSpeciesGCNorm$Recent))

meanBySpecies <- GC_VerVsNew %>%
	group_by(BinIndex, Species, Type) %>%
	summarise(meanGC = mean(GC_content, na.rm = TRUE)) %>%
	spread(Type, meanGC) %>%
	mutate(
		Presence = case_when(
			is.na(Recent) & !is.na(Ver) ~ "VerBin",
			!is.na(Recent) & !is.na(Ver) ~ "BothBin",
			!is.na(Recent) & is.na(Ver) ~ "RecentBin")
	) %>%
	filter(!is.na(Presence)) %>%
	gather(key = geneType, value = GC_norm, Ver, Recent, na.rm = TRUE) %>%
	mutate(geneType = factor(geneType, levels = c("Ver", "Recent"))) %>%
	mutate(
		Class = case_when(
			Presence == "VerBin" ~ "Vertical Only",
			Presence == "RecentBin" ~ "Recent Only",
			Presence == "BothBin" & geneType == "Ver" ~ "Vertical with Recent",
			Presence == "BothBin" & geneType == "Recent" ~ "Recent with Vertical"),
		Class = factor(Class, levels = c("Vertical Only", "Vertical with Recent", "Recent with Vertical", "Recent Only"))
	)

# Count the number of genes in each class
geneCounter_df	<- meanBySpecies %>% group_by(Class) %>% summarise(count = n(), meanPos = as.numeric(summary(GC_norm)[5]))

statComparisons		<- list(levels(meanBySpecies$Class)[1:2], levels(meanBySpecies$Class)[2:3], levels(meanBySpecies$Class)[3:4])
recentHGT_VerContext_boxplot4	<- ggplot(data = meanBySpecies, aes(x = Class, y = GC_norm, fill = geneType)) +
	scale_x_discrete(name = "Gene Class") +
	scale_y_continuous(name = "Normalised GC value") +
	geom_boxplot(color = axisCol, size = 0.7) +
	geom_label(data = geneCounter_df, aes(y = meanPos, label = count), color = textCol, fill = "white") +
	scale_fill_manual(values = c(dataTypeCols$Ver, dataTypeCols$Recent)) +
	stat_compare_means(comparisons = statComparisons, method = "wilcox.test", p.adjust = "bonferroni", color = axisCol, size = 0.8, label = "p.format") +
 	lightTheme +
 	theme(panel.grid.minor = element_blank())


quartz(width = 16, height = 10)
print(recentHGT_VerContext_boxplot)
quartz.save(file = file.path(GC_byBinFig_path, "recentGenesVerticalContext.pdf"), type = "pdf", dpi = 300)
invisible(dev.off())



# ## 767 of the 918 RecentOnly genes are in the 19 species with GI data
# recentOnly_data			<- subset(bySpecies_normGC_df, type == "RecentOnly")
# recentOnly_GISpecies	<- recentOnly_data[which(recentOnly_data$locusTag %in% byType_withGI_data$lHGT$`4`$locusTag),]
# GI_status	<- subset(byType_withGI_data$lHGT$`4`, locusTag %in% recentOnly_data$locusTag, select = In_GI, drop = TRUE)
# recentOnly_GISpecies$In_GI	<- GI_status
# nrow(subset(recentOnly_GISpecies, In_GI == TRUE))
# # 470 / 767 are in GIs

# ## 325 of the 401 RecentWithVer genes are in the 19 species with GI data
# recentWithVer_data		<- subset(bySpecies_normGC_df, type == "RecentWithVer")
# recentWithVer_GISpecies <- recentWithVer_data[which(recentWithVer_data$locusTag %in% byType_withGI_data$lHGT$`4`$locusTag),]
# GI_status	<- subset(byType_withGI_data$lHGT$`4`, locusTag %in% recentWithVer_data$locusTag, select = In_GI, drop = TRUE)
# recentWithVer_GISpecies$In_GI	<- GI_status
# nrow(subset(recentWithVer_GISpecies, In_GI == TRUE))
# # 110 / 325 are in GIs

# ## So, 61.2% of RecentOnly genes are in GIs - 33.9% of recentWithVer genes are in GIs, explaining the GC-content difference

# Write plot



# ------------------------------------------------------------------------------------- #
# Is there a correlation between the GC content of an lHGT gene and the vertical genes in the shared genomic bin?
# With recent transfers - are genes gained preferentially in more closely matching GC contexts?
# With older transfers - is there a general pattern of amelioration towards the vertical-gene GC context of it's environment?

GC_VerVsOldNew		<- bind_rows(list(byTypeSpeciesGCNorm$Ver, byTypeSpeciesGCNorm$Old, byTypeSpeciesGCNorm$Recent))

# Choose correlation method
corrMethod	<- "pearson"

verVsHGT_normGC_cor_list	<- lapply(c("Recent", "Old"), function(hgtType) {

	message(paste0("Working on ", hgtType, "..."))
	bySpecies_compare	<- mclapply(binomial_list, mc.cores = 14, function(speciesName) {

		# Iterate over bins
		byBinMeanGC_list	<- lapply(unique(GC_VerVsOldNew$BinIndex), function(genome_bin) {

			with(GC_VerVsOldNew, {
				byBin_HGT	<<- mean(GC_norm[which(BinIndex == genome_bin & Species == speciesName & Type == hgtType)])
				byBin_Ver	<<- mean(GC_norm[which(BinIndex == genome_bin & Species == speciesName & Type == "Ver")])
			})

			return(data.frame(BinIndex = genome_bin, Ver_GC = byBin_Ver, lHGT_GC = byBin_HGT, stringsAsFactors = FALSE))
		})
		byBinMeanGC_df	<- bind_rows(byBinMeanGC_list)

		bothPresent		<- subset(byBinMeanGC_df, !is.na(Ver_GC) & !is.na(lHGT_GC))
		if (nrow(bothPresent) == 0) bothPresent	<- data.frame(BinIndex = NA, Ver_GC = NA, lHGT_GC = NA, stringsAsFactors = FALSE)
		bothPresent$Species	<- speciesName

		return(bothPresent)
 
	})

	bySpecies_compare_df			<- bind_rows(bySpecies_compare)
	bySpecies_compare_df$comparison	<- hgtType

	# Calculate correlation
	calcCorr	<- rcorr(x = bySpecies_compare_df$Ver_GC, y = bySpecies_compare_df$lHGT_GC, type = corrMethod)
	corVal		<- signif(calcCorr$r[1,2], digits = 2)
	pVal		<- signif(calcCorr$P[1,2], digits = 2)

	# x + y positions
	x.pos		<- min(bySpecies_compare_df$Ver_GC, na.rm = TRUE)
	y.pos		<- mean(bySpecies_compare_df$lHGT_GC, na.rm = TRUE)

	# Cutoff for significance
	pval_lab	<- ifelse(pVal <= 0.05, paste0("~italic(p) == ", pVal), paste0("~NS"))

	# Values for plotting the correlation
	label		<- ifelse(identical(corrMethod, "pearson"), paste("italic(r) == ", corVal, pval_lab), paste("italic(rs) == ", corVal, pval_lab))
	pDf			<- data.frame(x.pos = x.pos, y.pos = y.pos, comparison = hgtType, corVal = corVal, pVal = pVal, plotLabel = label, stringsAsFactors = FALSE)

	return(list(allComparison_df = bySpecies_compare_df, statsComparison_df = pDf))
})

# Bind together the Ver vs Recent by Bin data and Ver vs Old by Bin data
verVsHGT_normGC_cor_df	<- bind_rows(lapply(verVsHGT_normGC_cor_list, function(element) return(element$allComparison_df)))

# Bind together the stats dataframes
verVsHGT_stats_df		<- bind_rows(lapply(verVsHGT_normGC_cor_list, function(element) return(element$statsComparison_df)))

# Set up colours
verVsHGT_normGC_cor_df$comparison	<- factor(verVsHGT_normGC_cor_df$comparison, levels = c("Old", "Recent"))

# Produce plot
verVsHGT_normGC_cor_plot	<- ggplot(data = verVsHGT_normGC_cor_df, aes(x = Ver_GC, y = lHGT_GC, color = comparison, fill = comparison)) +
	scale_x_continuous(name = "Vertical Gene normalised GC content", limits = c(-0.2, 0.1)) +
	scale_y_continuous(name = "HGT Gene normalised GC content") +
	geom_point(size = 0.6, alpha = 0.6) +
	geom_smooth(method = 'lm') +
	geom_abline(slope = 1, size = 1, intercept = 0, linetype = "longdash", color = axisCol) +
	geom_label(data = verVsHGT_stats_df, aes(x = x.pos, y = y.pos, label = plotLabel), vjust = 4, hjust = 1.2, parse = TRUE, size = 4, fill = "white", show.legend = FALSE) +
	scale_color_manual(values = c(dataTypeCols$Old, dataTypeCols$Recent), name = "Age of HGT") +
	scale_fill_manual(values = c(dataTypeCols$Old, dataTypeCols$Recent), name = "Age of HGT") +
	ggtitle("Comparing GC content of Vertical and HGT genes in the same genomic bins") +
	lightTheme

# Write plot
quartz(width = 16, height = 10)
print(verVsHGT_normGC_cor_plot)
quartz.save(file = file.path(GC_byBinFig_path, "byAgeHGTsVsVertical_GCcontent_byBin.pdf"), type = "pdf", dpi = 300)
invisible(dev.off())


# ------------------------------------------------------------------------------------- #

# verCrossSpecies	<- lapply(1:binNumber, function(genomeBin) {

# 	perBinGCNorm	<- subset(byTypeSpeciesGCNorm$Ver, BinIndex == genomeBin, select = GC_norm, drop = TRUE)
# 	av_GCNorm		<- mean(perBinGCNorm, na.rm = TRUE)
# 	out_df			<- data.frame(BinIndex = genomeBin, avGC = av_GCNorm, stringsAsFactors = FALSE)
# 	return(out_df)
# })

# verCrossSpecies_df		<- bind_rows(verCrossSpecies)

# verCrossSpecies_circ	<-  verCrossSpecies_df
# verCrossSpecies_circ$BinIndex	<- circular(x = (verCrossSpecies_circ$BinIndex / binNumber) * (2 * pi), zero = pi / 2, rotation = "clock", modulo = "2pi")

# lHGTDensity <- perTypeData$lHGT$'4'$circDensity
# allDensity	<- perTypeData$All$circDensity
# verDensity	<- perTypeData$Ver$'3'$circDensity

# numGenes	<- length(lHGTDensity$data)

# GC_color	<- wes_palette("BottleRocket2")[1]


# lines.circular(verCrossSpecies_circ$BinIndex, verCrossSpecies_df$avGC * 70, col = GC_color, plot.info = mainPlot, shrink = 10)
# lines.circular(verCrossSpecies_circ$BinIndex, rep(0, 200), col = alpha(GC_color, 0.5), plot.info = mainPlot, shrink = 10)















