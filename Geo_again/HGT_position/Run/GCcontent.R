#!/usr/bin/env Rscript

# Load master variables and HGT position functions
invisible(sapply(HGTPos.all, source, .GlobalEnv))

# Any libraries explicitly used in the script
require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load("RSQLite", "plyr", "dplyr","ggplot2", "ggtree", "wesanderson", "Biostrings", "ggpubr", "Hmisc")


# ------------------------------------------------------------------------------------- #
# Read in data

message("\nReading in data...", appendLF = FALSE)

# General position data
perTypeData			<- readRDS(file.path(positionData_path, "AG_perTypeData.rds"))

# Per-type COG data
perTypeCOG_data		<- readRDS(file.path(positionData_path, "AG_perTypeCOGData.rds"))

# Subgroup data
subgroupData		<- readRDS(file.path(positionData_path, "AG_subgroupData.rds"))

# Subdivision key
subDivisionKey_file	<- file.path(positionData_path, "subDivisionKeyData.rds")
if(!file.exists(subDivisionKey_file)) {
	stop("The subdivision data file \"", subDivisionKey_file, "\" is required.\nRun the subDivisionPrep.R script!")
} else {
	subDivisionKey_data	<- readRDS(file.path(positionData_path, "subDivisionKeyData.rds"))
	subDivisionKey_df	<- subDivisionKey_data$subDivisionKey_df
	subDivison_cols		<- subDivisionKey_data$subDivison_cols
}

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

# Quartz plotting options common to this script
quartz.options(canvas = "#333233", bg = "#333233", type = "pdf")


# ------------------------------------------------------------------------------------- #
# ---	Prepare phylogenies for plotting output - need tip order for plotting later	--- #
# ------------------------------------------------------------------------------------- #

# Get the binomial tree from the subgroup data set
AG_binomTime_tree	<- subgroupData$AG_binomTime_tree

# Plot with and without branch.lenght to get phylo and clado trees
binomTree_phylo <- ggtree(tr = AG_binomTime_tree, ladderize = FALSE, col =  "#D9D9D9", size = 1.5) +
	geom_tiplab(angle = -90, hjust = -0.05, col =  "#D9D9D9") +
	scale_x_reverse(limits = c(0.5, 0)) +
	coord_flip() +
	darkTheme +
	theme(
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank()
	)
binomTree_clado <- ggtree(tr = AG_binomTime_tree, branch.length = "none", ladderize = FALSE, col =  "#D9D9D9", size = 1.5) +
	geom_tiplab(angle = -90, hjust = -0.05, col =  "#D9D9D9") +
	scale_x_reverse(limits = c(20, 0)) +
	coord_flip() +
	darkTheme +
	theme(
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank()
	)

# Prepare tree as dataframe to extract tip order
treeAsDf	<- subset(fortify(AG_binomTime_tree, ladderize = FALSE), isTip)
tipOrder	<- with(treeAsDf, label[order(y, decreasing = T)])

# Write trees out to PDF for combination plots
quartz(width = 22, height = 12, file = file.path(GC_TreesFig_path, "BinomialTree_phylo.pdf"))
print(binomTree_phylo)
invisible(dev.off())

quartz(width = 16, height = 8, file = file.path(GC_TreesFig_path, "BinomialTree_clado.pdf"))
print(binomTree_clado)
invisible(dev.off())


# ------------------------------------------------------------------------------------- #
# ---	Prepare GC-content data structure - this is done over N (100) genomic bins	--- #
# ------------------------------------------------------------------------------------- #

message("Calculating GC content profiles for all dataTypes in parallel")
byTypebySpeciesGC_data	<- mclapply(dataTypes_withAge, mc.cores = 14, function(dataType) {

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

			# Get the protIDs and locusTags. NB <<- is assignment to global environment
			with(subsetData, {
				geneProtIDs	<<- protID[which(binomial == species & relGeneStart > bin_start & relGeneStart <= bin_end)]
				genelocTags	<<- locusTag[which(binomial == species & relGeneStart > bin_start & relGeneStart <= bin_end)]
			})

			# If the number of genes of this dataType in this Bin is 0 - return
			if (length(geneProtIDs) < geneCutOff) {
				return(data.frame(BinIndex = genome_bin, GC_content = NA, stringsAsFactors = FALSE))
			}

			# Get the corresponding nucleotide sequences from the Sqlite database
			nuclSeq_tbl		<- dbSendQuery(conn, 'SELECT NuclSeq FROM t1 WHERE protID = :protIDs AND locus = :locTags AND binomial = :species')
			dbBind(nuclSeq_tbl, param = list(protIDs = geneProtIDs, locTags = genelocTags, species = rep(species, length(geneProtIDs))))
			nuclSeq_df		<- dbFetch(nuclSeq_tbl)
			dbClearResult(nuclSeq_tbl)

			# Convert to DNAStringSet
			nuclSeq_stringSet	<- DNAStringSet(x = nuclSeq_df$NuclSeq)

			# For lHGT genes - calculate overall GC content
			GeneGC_cont		<- rowSums(subset(alphabetFrequency(nuclSeq_stringSet, as.prob = TRUE), select = 2:3, drop = FALSE))
			
			return(data.frame(BinIndex = genome_bin, GC_content = GeneGC_cont, stringsAsFactors = FALSE))
		})

		# Combine the per-Bin data
		byBinGC_df			<- bind_rows(byBinGC_list)
		byBinGC_df$Species	<- species

		return(byBinGC_df)
	})
	# Combine the per-Species data
	bySpeciesGC_df		<- bind_rows(bySpeciesGC_list)
	bySpeciesGC_df$Type	<- dataType

	return(bySpeciesGC_df)

})

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
	geom_violin(col = "#D9D9D9") +
	scale_fill_manual(values = perSpecies_cols, guide = FALSE) +
	darkTheme +
	theme(
		axis.text.x = element_blank()
	)

quartz(width = 16, height = 8, file = file.path(GC_SpeciesFig_path, "GC_allBySpecies.pdf"))
print(GC_all_plot)
invisible(dev.off())



# ---	Is there a difference in GC content between All lHGTs and Vertical genes?	--- #

GC_VervsHGT			<- bind_rows(list(byTypeSpeciesGCNorm$Ver, byTypeSpeciesGCNorm$lHGT))

# Factor variables and set up colours 
GC_VervsHGT$Type	<- factor(GC_VervsHGT$Type, levels = dataTypes_withAge)
GC_VervsHGT$Species	<- factor(GC_VervsHGT$Species, levels = rev(tipOrder))
verVsHGT_cols		<- wes_palette("Darjeeling1")[1:2]

# Prepare plot
GC_VervsHGT_plot	<- ggplot(data = GC_VervsHGT, aes(x = Type, y = GC_content, color = Type, fill = Type)) +
	geom_boxplot() +
	scale_y_continuous(limits = c(0.2, 0.7)) +
	facet_wrap(~Species, nrow = 1) +
	stat_compare_means(method = "t.test", col = "#D9D9D9", label = "p.signif", label.x.npc = "centre", label.y.npc = "bottom", size = 8) +
	scale_fill_manual(values = verVsHGT_cols) +
	darkTheme +
	theme(
		panel.spacing = unit(2, "lines"),
		axis.text.x = element_blank(),
		panel.grid.minor = element_blank(),
		panel.grid.major.x = element_blank(),
		axis.ticks.x = element_blank(),
		strip.background = element_rect(fill = "transparent", color = "#D9D9D9"),		
		strip.text = element_text(color =  "#D9D9D9")		
	)


quartz(width = 22, height = 6, file = file.path(GC_SpeciesFig_path, "GC_VervsHGTbySpecies.pdf"))
print(GC_VervsHGT_plot)
invisible(dev.off())


# ---	Is there a difference in GC content between OLD lHGTs and RECENT lHGTs?		--- #

# NB - we get some "Recent" values for species that don't have Recent ancestor branches
# This is due to secondary AG2AG transfers from other (recent) branches into these taxa

GC_OldvsNew			<- bind_rows(list(byTypeSpeciesGCNorm$Old, byTypeSpeciesGCNorm$Recent))

# Factor variables and set up colours
GC_OldvsNew$Type	<- factor(GC_OldvsNew$Type, levels = dataTypes_withAge)
GC_OldvsNew$Species	<- factor(GC_OldvsNew$Species, levels = rev(tipOrder))
oldVsNew_cols		<- wes_palette("Darjeeling1")[4:5]

# Here we count the number of observations for old and recent HGTs per species
bySpecies_count		<- lapply(binomial_list, function(speciesName) {

	byType_count	<- lapply(unique(GC_OldvsNew$Type), function(dataType) {
		byTypeData	<- subset(GC_OldvsNew, Species == speciesName & Type == dataType)
		max			<- max(byTypeData$GC_content, na.rm = TRUE)
		count		<- length(which(!is.na(byTypeData$GC_content)))
		return(data.frame(Type = dataType, Species = speciesName, Count = count, MaxVal = max, stringsAsFactors = FALSE))
	})

	return(bind_rows(byType_count))
})

bySpecies_count_df	<- bind_rows(bySpecies_count) 
bySpecies_count_df$Species	<- factor(bySpecies_count_df$Species, levels = rev(tipOrder))

# Plot
GC_OldvsNew_plot	<- ggplot(data = GC_OldvsNew, aes(x = Type, y = GC_content, color = Type, fill = Type)) +
	geom_boxplot() +
	scale_y_continuous(limits = c(0.2, 0.7)) +
	facet_wrap(~Species, nrow = 1) +
	geom_text(data = bySpecies_count_df, aes(y = MaxVal, label = Count), vjust = -1) +
	stat_compare_means(method = "t.test", col = "#D9D9D9", label = "p.signif", label.x.npc = "centre", label.y.npc = "bottom", size = 8) +
	scale_fill_manual(values = alpha(oldVsNew_cols, 0.5)) +
	scale_color_manual(values = oldVsNew_cols) +
	darkTheme + 
	theme(
		panel.spacing = unit(2, "lines"),
		axis.text.x = element_blank(),
		panel.grid.minor = element_blank(),
		panel.grid.major.x = element_blank(),
		axis.ticks.x = element_blank(),
		strip.background = element_rect(fill = "transparent", color = "#D9D9D9"),		
		strip.text = element_text(color =  "#D9D9D9")		
	)

# Write plot
quartz(width = 22, height = 6, file =  file.path(GC_SpeciesFig_path, "GC_OldvsNewbySpecies.pdf"))
print(GC_OldvsNew_plot)
invisible(dev.off())


# ---	Is there a difference in GC content between OLD lHGTs and Vertical genes?	--- #

GC_VervsOld			<- bind_rows(list(byTypeSpeciesGCNorm$Ver, byTypeSpeciesGCNorm$Old))

# Factor variables and set up colours
GC_VervsOld$Type	<- factor(GC_VervsOld$Type, levels = dataTypes_withAge)
GC_VervsOld$Species	<- factor(GC_VervsOld$Species, levels = rev(tipOrder))
VervsOld_cols		<- wes_palette("Darjeeling1")[c(1,4)]

bySpecies_count		<- lapply(binomial_list, function(speciesName) {

	byType_count	<- lapply(unique(GC_VervsOld$Type), function(dataType) {
		byTypeData	<- subset(GC_VervsOld, Species == speciesName & Type == dataType)
		max			<- max(byTypeData$GC_content, na.rm = TRUE)
		count		<- length(which(!is.na(byTypeData$GC_content)))
		return(data.frame(Type = dataType, Species = speciesName, Count = count, MaxVal = max, stringsAsFactors = FALSE))
	})

	return(bind_rows(byType_count))
})
bySpecies_count_df	<- bind_rows(bySpecies_count)
bySpecies_count_df$Species	<- factor(bySpecies_count_df$Species, levels = rev(tipOrder))


GC_VervsOld_plot	<- ggplot(data = GC_VervsOld, aes(x = Type, y = GC_content, color = Type, fill = Type)) +
	geom_boxplot() +
	scale_y_continuous(limits = c(0.2, 0.7)) +
	facet_wrap(~Species) +
	geom_text(data = bySpecies_count_df, aes(y = MaxVal, label = Count), vjust = -1) +
	stat_compare_means(method = "t.test", col = "#D9D9D9", label = "p.signif", label.x.npc = "centre", label.y.npc = "bottom", size = 8) +
	scale_fill_manual(values = alpha(VervsOld_cols, 0.5)) +
	scale_color_manual(values = VervsOld_cols) +
	darkTheme + 
	theme(
		panel.spacing = unit(2, "lines"),
		axis.text.x = element_blank(),
		panel.grid.minor = element_blank(),
		panel.grid.major.x = element_blank(),
		axis.ticks.x = element_blank(),
		strip.background = element_rect(fill = "transparent", color = "#D9D9D9"),		
		strip.text = element_text(color =  "#D9D9D9")		
	)

# Write plot
quartz(width = 14, height = 18, file =  file.path(GC_SpeciesFig_path, "GC_VervsOldbySpecies.pdf"))
print(GC_VervsOld_plot)
invisible(dev.off())






# ------------------------------------------------------------------------------------- #
# ---	Examine GC-content across the genome by looking at perBin GC differences	--- #
# ------------------------------------------------------------------------------------- #


# ---		Convert the subdivision structure to a linear form for plotting			--- #

xValCols		<- subset(subDivisionKey_df, select = c(xmin, xmax))
xValLineCols	<- 100 * (xValCols / (2 * pi))
subDivisionKeyLinear_df	<- cbind(xValLineCols, subset(subDivisionKey_data$subDivisionKey_df, select = -c(xmin, xmax)))


# ---	Compare all-Gene GC-content change per-Species across the genome Ori-to-Ori	--- #

# Get the relevant data + factor
GC_allGenes				<- byTypeSpeciesGCNorm$All
GC_allGenes$BinIndex	<- (GC_allGenes$BinIndex * 100) / binNumber
GC_allGenes$Type	<- factor(GC_allGenes$Type, levels = dataTypes_withAge)

# Reassign ymin/ymax for subdivision bars for the allGenes set
allGenes_subdivs		<- subDivisionKeyLinear_df
allGenes_subdivs$ymin	<- round_any(min(GC_allGenes$GC_content, na.rm = TRUE), accuracy = 0.1, f = ceiling)
allGenes_subdivs$ymax	<- round_any(max(GC_allGenes$GC_content, na.rm = TRUE), accuracy = 0.1, f = floor)

# Color palette for allGenes
allGenes_col		<- "#D9D9D9"

# Plot. This plot is not very informative, as there is a lot of intraGenome fluctuations and the two GC-distinct groups complicate viewing
GCbyBin_all_plot	<- ggplot(GC_allGenes, aes(x = BinIndex, y = GC_content, color = Type)) +
	scale_x_continuous(name = "Percentage along genome", limits = c(0, 100), breaks = seq(0, 100, 10), minor_breaks = NULL) +
	scale_y_continuous(name = "GC content") +
	geom_rect(data =  allGenes_subdivs, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = alpha(unlist(lapply(subDivison_cols, rep, times = 2)), 0.2), inherit.aes = FALSE) +
	geom_smooth(
		data = GC_allGenes,
		aes(group = Species),
		method = "loess",
		span = 0.1,
		col = alpha(allGenes_col, 0.8),
		size = 0.2,
		se = FALSE) +
	# Overall smooth line (cross-species)
	scale_color_manual(values = allGenes_col, guide = guide_legend(title = "Gene Type")) +
	scale_fill_manual(values = alpha(allGenes_col, 0.7), guide = FALSE) +
	darkTheme

# Write plot
quartz(width = 12, height = 8, dpi = 300, type = "png", file = file.path(GC_byBinFig_path, "perSpeciesAll_raw.png"))
print(GCbyBin_all_plot)
invisible(dev.off())


# ---		Compare normalised GC-content change per-Species across the genome 		--- #

# Set the new subdivision dimensions
allGenesNorm_subdivs		<- subDivisionKeyLinear_df
allGenesNorm_subdivs$ymin	<- -0.15
allGenesNorm_subdivs$ymax	<- 0.15

# Plot with normalised per-Bin GC values. Have to use GAM smoothhing method (loess memory fail), k=25 looks like happy medium.
GCbyBin_allNorm_plot	<- ggplot(GC_allGenes, aes(x = BinIndex, y = GC_norm, color = Type)) +
	scale_x_continuous(name = "Percentage along genome", limits = c(0, 100), breaks = seq(0, 100, 10), minor_breaks = NULL) +
	scale_y_continuous(name = "GC content") +
	geom_rect(data =  allGenesNorm_subdivs, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = alpha(unlist(lapply(subDivison_cols, rep, times = 2)), 0.2), inherit.aes = FALSE) +
	geom_smooth(
		data = GC_allGenes,
		aes(group = Species),
		method = "loess",
		span = 0.1,
		col = alpha(allGenes_col, 0.6),
		size = 0.2,
		se = FALSE) +
	# Overall smooth line (cross-species)
	geom_smooth(method = "gam", formula = y ~ s(x, k = 25), size = 1.5, se = TRUE, span = 0.2, color = "red") +
	scale_color_manual(values = allGenes_col, guide = guide_legend(title = "Gene Type")) +
	scale_fill_manual(values = alpha(allGenes_col, 0.7), guide = FALSE) +
	darkTheme

# Write plot
quartz(width = 12, height = 8, dpi = 300, type = "png", file = file.path(GC_byBinFig_path, "perSpeciesAll_normalised.png"))
print(GCbyBin_allNorm_plot)
invisible(dev.off())


# ---		Compare normalised GC-content between Ver and lHGT genes per-Species	--- #

# Get the relevant data
GC_VerVsHGT				<- bind_rows(list(byTypeSpeciesGCNorm$Ver, byTypeSpeciesGCNorm$lHGT))
GC_VerVsHGT$BinIndex	<- (GC_VerVsHGT$BinIndex * 100) / binNumber

# Factor variables + color
GC_VerVsHGT$Type	<- factor(GC_VerVsHGT$Type, levels = dataTypes_withAge)
GC_VerVsHGT$Species	<- factor(GC_VerVsHGT$Species, levels = rev(tipOrder))
verVsHGT_cols		<- wes_palette("Darjeeling1")[1:2]

# Plot with normalised per-Bin GC values. Have to use GAM smoothhing method (loess memory fail), k=20 looks most like loess.
GCbyBin_VerHGT_norm_plot	<- ggplot(GC_VerVsHGT, aes(x = BinIndex, y = GC_norm, color = Type), fill = Type) +
	scale_x_continuous(name = "Percebtage along genome", limits = c(0, 100), breaks = seq(0, 100, 10), minor_breaks = NULL) +
	scale_y_continuous(name = "Normalised GC content") +
	geom_rect(data =  allGenesNorm_subdivs, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = alpha(unlist(lapply(subDivison_cols, rep, times = 2)), 0.2), inherit.aes = FALSE) +
	geom_smooth(
		data = subset(GC_VerVsHGT, Type == "Ver"),
		aes(group = Species),
		method = "loess",
		span = 0.1,
		col = alpha(verVsHGT_cols[2], 0.4),
		size = 0.2,
		se = FALSE) +
	geom_smooth(
		data = subset(GC_VerVsHGT, Type == "lHGT"),
		aes(group = Species),
		method = "loess",
		span = 0.1,
		col = alpha(verVsHGT_cols[1], 0.4),
		size = 0.2,
		se = FALSE) +
	# Overall smooth line (cross-species)
	geom_smooth(aes(color = Type, fill = Type), method = "gam", formula = y ~ s(x, k = 20), size = 1.5, se = TRUE, span = 0.1) +
	scale_color_manual(values = verVsHGT_cols, guide = guide_legend(title = "Gene Type")) +
	scale_fill_manual(values = alpha(verVsHGT_cols, 0.7), guide = FALSE) +
	darkTheme +
	theme(
		legend.justification = c(1, 1),
		legend.position = c(0.9, 0.9),
		legend.background = element_rect(fill = "#333233", color = "#D9D9D9")
	)

# Write plot
quartz(width = 18, height = 10, dpi = 300, type = "png", file = file.path(GC_byBinFig_path, "perSpeciesVerHGT_normalised.png"))
print(GCbyBin_VerHGT_norm_plot)
invisible(dev.off())




# ---	Compare Vertical v lHGT GC-content per-Species across the genome Ori-to-Ori	--- #

# Reassign ymin/ymax for the subdivision bars for the GC content scale. NB loess warnings are unimportant
verVsHGT_subdivs	<- lapply(binomial_list, function(speciesName) {
	bySpecies		<- subset(GC_VerVsHGT, Species == speciesName)
	byType_list		<- lapply(list("Ver", "lHGT"), function(dataType) {
		byType		<- subset(bySpecies, Type == dataType)
		byTypeLoess	<- loess.smooth(byType$BinIndex, byType$GC_content, span = 0.1, degree = 1)
		return(data.frame(min = min(byTypeLoess$y), max = max(byTypeLoess$y)))
	})
	byType_df		<- bind_rows(byType_list)
	
	subDivisionKeyLinear_df$ymin	<- round_any(min(byType_df$min), accuracy = 0.1, f = floor)
	subDivisionKeyLinear_df$ymax	<- round_any(max(byType_df$max), accuracy = 0.1, f = ceiling)
	subDivisionKeyLinear_df$Species	<- speciesName
	return(subDivisionKeyLinear_df)
})
verVsHGT_subdivs_df	<- bind_rows(verVsHGT_subdivs)

# Color palette for verVslHGT
verVsHGT_cols		<- wes_palette("Darjeeling1")[1:2]

# Plot
GCbySpecies_VerVsHGT_plot	<- ggplot(GC_VerVsHGT, aes(x = BinIndex, y = GC_content, color = Type)) +
	scale_x_continuous(name = "Percentage along genome", limits = c(0, 100), breaks = NULL, minor_breaks = NULL) +
	scale_y_continuous(name = "GC content") +
	geom_rect(data =  verVsHGT_subdivs_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, group = Species), fill = rep(alpha(unlist(lapply(subDivison_cols, rep, times = 2)), 0.2), 25), inherit.aes = FALSE) +
	geom_smooth(
		data = subset(GC_VerVsHGT, Type == "lHGT"),
		aes(group = Species),
		method = "loess",
		span = 0.2,
		color = alpha(verVsHGT_cols[1], 1),
		size = 0.2,
		se = FALSE) +
	geom_smooth(
		data = subset(GC_VerVsHGT, Type == "Ver"),
		aes(group = Species),
		method = "loess",
		span = 0.2,
		color = alpha(verVsHGT_cols[2], 1),
		size = 0.2,
		se = FALSE) +
	facet_wrap(~Species, scales = "free") +
	scale_color_manual(values = verVsHGT_cols, guide = FALSE) +
	scale_fill_manual(values = alpha(verVsHGT_cols, 0.7), guide = FALSE) +
	darkTheme +
	theme(
		panel.spacing = unit(2, "lines"),
		panel.grid.minor = element_blank(),
		strip.background = element_rect(fill = "transparent", color = "#D9D9D9"),		
		strip.text = element_text(color =  "#D9D9D9")
	)

# Write plot
quartz(width = 24, height = 14, dpi = 300, type = "png", file = file.path(GC_byBinFig_path, "bySpeciesVerHGT_facet.png"))
print(GCbySpecies_VerVsHGT_plot)
invisible(dev.off())




# ---	Correlate Vertical & lHGT GC-content by Species over the genome Ori-to-Ori	--- #

# Get the relevant data - don't use the one above (normalised to 100% x-axis)
GC_VerVsHGT				<- bind_rows(list(byTypeSpeciesGCNorm$Ver, byTypeSpeciesGCNorm$lHGT))

# Choose correlation method
corrMethod	<- "pearson"

# For each species, for each bin, calculate the mean GC content for the present Vertical and lHGT genes
# Correlate these genome wide - to see whether there is an association
bySpecies_meanGC_list	<- mclapply(binomial_list, mc.cores = 14, function(speciesName) {

	byBinMeanGC_list	<- lapply(1:binNumber, function(genome_bin) {

		with(GC_VerVsHGT, {
			byBin_HGT	<<- mean(GC_content[which(BinIndex == genome_bin & Species == speciesName & Type == "lHGT")])
			byBin_Ver	<<- mean(GC_content[which(BinIndex == genome_bin & Species == speciesName & Type == "Ver")])
		})

		return(data.frame(BinIndex = genome_bin, Ver_GC = byBin_Ver, lHGT_GC = byBin_HGT, stringsAsFactors = FALSE))
	})
	byBinMeanGC_df	<- bind_rows(byBinMeanGC_list)
	byBinMeanGC_df$Species	<- speciesName

	# Correlation over bins that contain non-missing GC content
	nonMissing	<- byBinMeanGC_df[which(rowSums(is.na(byBinMeanGC_df)) == 0),]
	calcCorr	<- rcorr(x = nonMissing[,2], y = nonMissing[,3], type = corrMethod)
	corVal		<- signif(calcCorr$r[1,2], digits = 2)
	pVal		<- signif(calcCorr$P[1,2], digits = 2)

	# Cutoff for significance
	pval_lab	<- ifelse(pVal <= 0.05, paste0("~italic(p) == ", pVal), paste0("~NS"))

	# Values for plotting the correlation
	label		<- ifelse(identical(corrMethod, "pearson"), paste("italic(r) == ", corVal, pval_lab), paste("italic(rs) == ", corVal, pval_lab))
	pos.x		<- mean(range(nonMissing$Ver_GC))
	pos.y		<- max(nonMissing$lHGT_GC) - (0.05 * max(nonMissing$lHGT_GC))
	pDf			<- data.frame(x = pos.x, y = pos.y, nBinsIncl = nrow(nonMissing), corVal = corVal, pVal = pVal, Species = speciesName, plotLabel = label, stringsAsFactors = FALSE)

	return(list(GCobs = byBinMeanGC_df, cor = pDf))
})

# Prepare the correlation and label data
bySpecies_meanGC_df	<- bind_rows(lapply(bySpecies_meanGC_list, function(l) l$GCobs))
bySpecies_PCorr_df	<- bind_rows(lapply(bySpecies_meanGC_list, function(l) l$cor))

# Factor the data and labels for the correlation plot
bySpecies_meanGC_df$Species	<- factor(bySpecies_meanGC_df$Species, levels = rev(tipOrder))
bySpecies_PCorr_df$Species	<- factor(bySpecies_PCorr_df$Species, levels = rev(tipOrder))

# Plot the correlation plot
VervHGT_GCCorr_plot	<- ggplot(data = bySpecies_meanGC_df, aes(x = Ver_GC, y = lHGT_GC, group = Species)) + 
	scale_y_continuous(name = "lHGT GC Content per-Bin") +
	scale_x_continuous(name = "Vertical GC Content per-Bin", limits = c(0.3, 0.6)) +
	geom_point(size = 1, color = "#D9D9D9") +
	geom_smooth(method = "lm", color = "red") +
	facet_wrap(~Species) +
	# With scales free, use the pos.x and pos.y values
	geom_text(data = bySpecies_PCorr_df, aes(x = 0.45, y = 0.55, label = plotLabel), parse = TRUE, size = 4, color = "red") +
	darkTheme +
	theme(
		panel.spacing = unit(2, "lines"),
		panel.grid.minor = element_blank(),
		strip.background = element_rect(fill = "transparent", color = "#D9D9D9"),		
		strip.text = element_text(color =  "#D9D9D9")
	)

# Which species have significantly correlated Ver and HGT content across genomic bins?
VerHGT_corrSpecies	<- bySpecies_PCorr_df[which(bySpecies_PCorr_df$pVal < 0.01),]

# Write plot
quartz(width = 24, height = 14, dpi = 300, type = "png", file = file.path(GC_byBinFig_path, "correlateBySpecies_VerHGT_facet.png"))
print(VervHGT_GCCorr_plot)
invisible(dev.off())





# ------------------------------------------------------------------------------------- #
# Does the genomic context of HGT gain location affect either the vertical genes or HGT genes gained?
# Compare the GC content of recent HGT gains across genomes (normalised GC), in bins where there are vertical
# genes and bins where there are no vertical genes

# Need Vertical and recent HGT data
GC_VerVsNew		<- bind_rows(list(byTypeSpeciesGCNorm$Ver, byTypeSpeciesGCNorm$Recent))
valueTypes		<- c("VerOnly", "VerWithRecent", "RecentWithVer", "RecentOnly")

# Sort data depending on the content of genes in the bin
bySpecies_normGC_list	<- mclapply(binomial_list, mc.cores = 14, function(speciesName) {

	byBinMeanGC_list	<- lapply(1:binNumber, function(genome_bin) {

		with(GC_VerVsNew, {
			byBin_HGT	<<- mean(GC_norm[which(BinIndex == genome_bin & Species == speciesName & Type == "Recent")])
			byBin_Ver	<<- mean(GC_norm[which(BinIndex == genome_bin & Species == speciesName & Type == "Ver")])
		})

		return(data.frame(BinIndex = genome_bin, Ver_GC = byBin_Ver, lHGT_GC = byBin_HGT, stringsAsFactors = FALSE))
	})
	byBinMeanGC_df	<- bind_rows(byBinMeanGC_list)

	byTypeGC_cont	<- lapply(valueTypes, function(subsetType) {

		if (identical(subsetType, "VerOnly")) {
			values	<- subset(byBinMeanGC_df, !is.na(Ver_GC) & is.na(lHGT_GC), select = Ver_GC, drop = TRUE)
		} else if (identical(subsetType, "VerWithRecent")) {
			values	<- subset(byBinMeanGC_df, !is.na(Ver_GC) & !is.na(lHGT_GC), select = Ver_GC, drop = TRUE)
		} else if (identical(subsetType, "RecentWithVer")) {
			values	<- subset(byBinMeanGC_df, !is.na(Ver_GC) & !is.na(lHGT_GC), select = lHGT_GC, drop = TRUE)
		} else if (identical(subsetType, "RecentOnly")) {	
			values	<- subset(byBinMeanGC_df, is.na(Ver_GC) & !is.na(lHGT_GC), select = lHGT_GC, drop = TRUE)
		}

		if (length(values) == 0) values	<- NA

		return(data.frame(GC_value = values, type = subsetType, stringsAsFactors = FALSE))
	})
	byTypeGC_cont	<- bind_rows(byTypeGC_cont)
	byTypeGC_cont$Species	<- speciesName

	return(byTypeGC_cont)
})
bySpecies_normGC_df	<- bind_rows(bySpecies_normGC_list)


# Factor and stat comparisons
bySpecies_normGC_df$type	<- factor(bySpecies_normGC_df$type, levels = valueTypes)
statComparisons				<- list(valueTypes[1:2], valueTypes[2:3], valueTypes[3:4])

# Counter for the number of bins in each category
binCounter_df	<- as.data.frame(table(bySpecies_normGC_df$type))
names(binCounter_df)	<- c("type", "Count")
binCounter_df$yval		<- round_any(max(bySpecies_normGC_df$GC_value, na.rm = T), accuracy = 0.01, f = floor)
binCounter_df$label		<- paste0("Bins = ", binCounter_df$Freq)

# Plot and plot colours
plotCols	<- wes_palette("Rushmore1")[2:5]
recentHGT_VerContext_boxplot	<- ggplot(data = bySpecies_normGC_df, aes(x = type, y = GC_value, fill = type)) +
	scale_y_continuous(name = "Normalised GC value") +
	scale_x_discrete(name = "Gene type", labels = str_wrap(c("Vertical Only", "Vertical with recent HGTs", "Recent HGTs with Vertical", "Recent HGTs Only"), width = 20)) +
	geom_boxplot(color = "#D9D9D9") +
	geom_label(data = binCounter_df, aes(y = yval, label = Count), vjust = 1, color = "#D9D9D9", fill = "#333233") +
	scale_fill_manual(values = plotCols) +
	stat_compare_means(mapping = aes(x = type, y = GC_value), comparisons = statComparisons, method = "wilcox.test", p.adjust = "bonferroni", color = "#D9D9D9", size = 0.5, label = "p.format") +
	darkTheme


# Write plot
quartz(width = 16, height = 10, dpi = 300, type = "png", file = file.path(GC_byBinFig_path, "recentGenesVerticalContext.png"))
print(recentHGT_VerContext_boxplot)
invisible(dev.off())


# ------------------------------------------------------------------------------------- #
# Is there a correlation between the GC content of an lHGT gene and the vertical genes in the shared genomic bin?
# With recent transfers - are genes gained preferentially in more closely matching GC contexts?
# With older transfers - is there a general pattern of amelioration towards the vertical-gene GC context of it's environment?

GC_VerVsHGT		<- bind_rows(list(byTypeSpeciesGCNorm$Ver, byTypeSpeciesGCNorm$Old, byTypeSpeciesGCNorm$Recent))

# Choose correlation method
corrMethod	<- "pearson"

verVsHGT_normGC_cor_list	<- lapply(c("Recent", "Old"), function(hgtType) {

	message(paste0("Working on ", hgtType, "..."))
	bySpecies_compare	<- mclapply(binomial_list, mc.cores = 14, function(speciesName) {

		# Iterate over bins
		byBinMeanGC_list	<- lapply(1:binNumber, function(genome_bin) {

			with(GC_VerVsHGT, {
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
comparisonCols			<- wes_palette("GrandBudapest1")[1:2]

# Produce plot
verVsHGT_normGC_cor_plot	<- ggplot(data = verVsHGT_normGC_cor_df, aes(x = Ver_GC, y = lHGT_GC, color = comparison, fill = comparison)) +
	scale_x_continuous(name = "Vertical Gene normalised GC content", limits = c(-0.2, 0.1)) +
	scale_y_continuous(name = "HGT Gene normalised GC content") +
	geom_point(size = 0.5, alpha = 0.4) +
	geom_smooth(method = 'lm') +
	geom_abline(slope = 1, size = 1, intercept = 0, linetype = "longdash", color = wes_palette("IsleofDogs1")[5]) +
	geom_label(data = verVsHGT_stats_df, aes(x = x.pos, y = y.pos, label = plotLabel), vjust = 4, hjust = 1.2, parse = TRUE, size = 4, fill =  "#333233", show.legend = FALSE) +
	scale_color_manual(values = comparisonCols, name = "Age of HGT") +
	scale_fill_manual(values = comparisonCols, name = "Age of HGT") +
	ggtitle("Comparing GC content of Vertical and HGT genes in the same genomic bins") +
	darkTheme

# Write plot
quartz(width = 16, height = 10, dpi = 300, type = "png", file = file.path(GC_byBinFig_path, "byAgeHGTsVsVertical_GCcontent_byBin.png"))
print(verVsHGT_normGC_cor_plot)
invisible(dev.off())


# ------------------------------------------------------------------------------------- #

verCrossSpecies	<- lapply(1:binNumber, function(genomeBin) {

	perBinGCNorm	<- subset(byTypeSpeciesGCNorm$Ver, BinIndex == genomeBin, select = GC_norm, drop = TRUE)
	av_GCNorm		<- mean(perBinGCNorm, na.rm = TRUE)
	out_df			<- data.frame(BinIndex = genomeBin, avGC = av_GCNorm, stringsAsFactors = FALSE)
	return(out_df)
})

verCrossSpecies_df		<- bind_rows(verCrossSpecies)

verCrossSpecies_circ	<-  verCrossSpecies_df
verCrossSpecies_circ$BinIndex	<- circular(x = (verCrossSpecies_circ$BinIndex / binNumber) * (2 * pi), zero = pi / 2, rotation = "clock", modulo = "2pi")

lHGTDensity <- perTypeData$lHGT$'4'$circDensity
allDensity	<- perTypeData$All$circDensity
verDensity	<- perTypeData$Ver$'3'$circDensity

numGenes	<- length(lHGTDensity$data)

GC_color	<- wes_palette("BottleRocket2")[1]


lines.circular(verCrossSpecies_circ$BinIndex, verCrossSpecies_df$avGC * 70, col = GC_color, plot.info = mainPlot, shrink = 10)
lines.circular(verCrossSpecies_circ$BinIndex, rep(0, 200), col = alpha(GC_color, 0.5), plot.info = mainPlot, shrink = 10)















