#!/usr/bin/env Rscript

# Load master variables and HGT position functions
invisible(sapply(HGTPos.all, source, .GlobalEnv))

# Any libraries explicitly used in the script
require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load("dplyr", "reshape2", "grid", "ggplot2", "ggrepel", "ggdendro", "ggpubr", "wesanderson")

# ------------------------------------------------------------------------------------- #
message("\nReading in data...", appendLF = FALSE)

# General position data
perTypeData			<- readRDS(file.path(positionData_path, "AG_perTypeData.rds"))

# Subgroup data
dnaA_pos_data		<- readRDS(file.path(positionData_path, "bySpecies_dnaA_data.rds"))

message("\rReading in data... done\n")
# ------------------------------------------------------------------------------------- #
# Open All_prot database
dbConn			<- dbConnect(RSQLite::SQLite(), allProtDB_path)

binomial_list		<- unique(perTypeData$All$allPosData$binomial)

# ------------------------------------------------------------------------------------- #
# Path to sporulation analysis
sporulation_path	<- file.path(master_path, "Sporulation")

# Read in the sigF regulon genes (Wang JMB 2006)
bsub_sigF_regulonGenes_file	<- file.path(sporulation_path, "bsubPY79_sigF_regulonGenes.tsv")
bsub_sigF_regulonGenes_data	<- read.table(bsub_sigF_regulonGenes_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
names(bsub_sigF_regulonGenes_data)		<- c("bSub_sigF_genes", "bSub_sigF_protIDs")

# B subtilis 168 taxid
bsub_168_taxid	<- 224308

# Find the protIDs for these genes in the host genome (all protIDs found for B sub 168; taxid = 224308).
bsub_sigF_regulonGenes_data$sigF_protIDs_wild	<- paste0(bsub_sigF_regulonGenes_data$bSub_sigF_protIDs, "*")

# Find the full prot IDs for these proteins in our database
bsub_sigF_protID_tbl	<- dbSendQuery(dbConn, 'SELECT protID, gene_start, OrthGroup FROM t1 WHERE taxid == :bsub_taxid AND protID GLOB :productsWild')
dbBind(bsub_sigF_protID_tbl, param = list(bsub_taxid = rep(bsub_168_taxid, nrow(bsub_sigF_regulonGenes_data)), productsWild = bsub_sigF_regulonGenes_data$sigF_protIDs_wild))
bsub_sigF_protID_df	<- dbFetch(bsub_sigF_protID_tbl)
dbClearResult(bsub_sigF_protID_tbl)

bsub_sigF_full_df	<- cbind(bsub_sigF_regulonGenes_data[,1:2], bsub_sigF_protID_df)

# Get the genome length for B subtilis 168
bsub_168_genomeLen_tmp	<- dbSendQuery(dbConn, 'SELECT genome_l FROM t1 WHERE taxid == 224308 LIMIT 1')
bsub_168_genomeLength	<- dbFetch(bsub_168_genomeLen_tmp)
dbClearResult(bsub_168_genomeLen_tmp)
bsub_168_genomeLength	<- bsub_168_genomeLength[1,1]

# Get the dnaA data for B subtilis 168
bsub_168_dnaA_data		<- subset(dnaA_pos_data, taxid == bsub_168_taxid)

# Get the relative genome positions
bsub_sigF_relStart_list	<- lapply(bsub_sigF_full_df$gene_start, genomeRelativePosition, oriStart = bsub_168_dnaA_data$oriStart, oriEnd = bsub_168_dnaA_data$oriEnd, oriStrand = bsub_168_dnaA_data$oriStrand, genomeLength = bsub_168_genomeLength)
bsub_sigF_full_df$relStart		<- unlist(bsub_sigF_relStart_list)
bsub_sigF_full_df$toOriStart	<- ifelse(bsub_sigF_full_df$relStart > 0.5, 1 - bsub_sigF_full_df$relStart, bsub_sigF_full_df$relStart)
bsub_sigF_full_df	<- bsub_sigF_full_df[order(-bsub_sigF_full_df$toOriStart),]

# Sample plot for G sub (compare to Wang JMB 2006 Fig 1)
bsub_sigF_genePosition_plot		<- ggplot(bsub_sigF_full_df, aes(x = relStart, xend = relStart, y = 1, yend = 0.8, label = bSub_sigF_genes)) +
	scale_y_continuous(limits = c(0, 1.2), breaks = 1.0) +
	scale_x_continuous(limits = c(0, 1), breaks = seq(0, 0.9, by = 0.5), labels = c("Origin", "Terminus")) +
	geom_segment(color = alpha(wes_palette("Zissou1")[5], 0.6)) +
	coord_polar() +
	geom_text_repel(size = 4, min.segment.length = 2, point.padding = 0.1) + 
	theme_bw()

# ------------------------------------------------------------------------------------- #
# Of the 48 genes from Wang, 43 are represented in Geobacillus orthologous groups
# The 43 genes belong to seperate orthologous families in AG

# sigF_gene_missingAG	<- bsub_sigF_full_df[is.na(bsub_sigF_full_df$OrthGroup),]
# sigF_gene_presentAG	<- bsub_sigF_full_df[!is.na(bsub_sigF_full_df$OrthGroup),]

# In 32 of the 43 orthologous groups, all taxa were represented
sigF_inAG_data <- lapply(1:nrow(bsub_sigF_full_df), function(sigF_index) {
	sigF_gene	<- bsub_sigF_full_df[sigF_index,]

	sigF_orthGroup	<- sigF_gene$OrthGroup
	sigF_geneName	<- sigF_gene$bSub_sigF_genes
	message(sigF_orthGroup)

	if (is.na(sigF_orthGroup)) {
		emptyRow			<- subset(perTypeData$All$allPosData[1,], select = -c(CircStart, CircEnd))
		emptyRow[1,]		<- NA
		emptyRow$geneName	<- sigF_geneName
		emptyRow$inAllAG	<- FALSE
		return(list(perOrth_data = emptyRow, uniqueTaxa = 0))
	}

	AG_inOrth	<- subset(perTypeData$All$allPosData, orthGroup == sigF_orthGroup & !taxid %in% outlierTaxid, select = -c(CircStart, CircEnd))
	AG_inOrth$geneName	<- sigF_geneName

	uniqueTaxa	<- length(unique(AG_inOrth$binomial))
	if (uniqueTaxa == 23) {
		presentInAll	<- TRUE
	} else {
		presentInAll	<- FALSE
	}

	AG_inOrth$inAllAG	<- presentInAll

	return(list(perOrth_data = AG_inOrth, uniqueTaxa = uniqueTaxa))
})

sigF_inAG_df	<- bind_rows(lapply(sigF_inAG_data, function(x) return(x$perOrth_data)))
sigF_inAG_count	<- unlist(lapply(sigF_inAG_data, function(x) return(x$uniqueTaxa)))
numGenes		<- length(unique(sigF_inAG_df$geneName))

siG_allInAG_df	<- subset(sigF_inAG_df, inAllAG == TRUE)
numInAllAG		<- length(unique(siG_allInAG_df$geneName))

# Just all genes where an ortholog is present in at least one AG species
sigF_inAG_plot_list	<- lapply(1:numGenes, function(geneIndex) {
	gene		<- unique(sigF_inAG_df$geneName)[geneIndex]
	geneData	<- subset(sigF_inAG_df, geneName == gene)

	geneData$xmin	<- geneIndex - 0.9
	geneData$xmax	<- geneIndex
	geneData$toOriStart	<- ifelse(geneData$relGeneStart > 0.5, 1 - geneData$relGeneStart, geneData$relGeneStart)

	bsubRelStart	<- subset(bsub_sigF_full_df, bSub_sigF_genes == gene, select = relStart, drop = TRUE)
	geneData$bsubRelStart	<- bsubRelStart
	geneData$bsubToOriStart	<- ifelse(bsubRelStart > 0.5, 1 - bsubRelStart, bsubRelStart)

	if (bsubRelStart <= 0.5) replichore <- "Right" else replichore <- "Left"
	geneData$bsubReplichore	<- replichore

	return(geneData)
})
sigF_inAG_plot_df		<- bind_rows(sigF_inAG_plot_list)


bSubOnly_pos_list	<- lapply(unique(sigF_inAG_plot_df$geneName), function(gene) {
	perGene_data	<- subset(sigF_inAG_plot_df, geneName == gene)
	singleEntry		<- perGene_data[1,]
	singleEntry[,c(1:12, 17)]	<- NA
	return(singleEntry)
})
bSubOnly_pos_df		<- bind_rows(bSubOnly_pos_list)

# # Just for those genes in which each AG species has an ortholog
# sigF_allInAG_plot_list	<- lapply(1:numInAllAG, function(geneIndex) {
# 	gene		<- unique(siG_allInAG_df$geneName)[geneIndex]
# 	geneData	<- subset(siG_allInAG_df, geneName == gene)

# 	geneData$xmin	<- geneIndex - 0.9
# 	geneData$xmax	<- geneIndex
# 	geneData$toOriStart	<- ifelse(geneData$relGeneStart > 0.5, 1 - geneData$relGeneStart, geneData$relGeneStart)

# 	bsubRelStart	<- subset(bsub_sigF_full_df, bSub_sigF_genes == gene, select = relStart, drop = TRUE)
# 	geneData$bsubRelStart	<- bsubRelStart
# 	geneData$bsubToOriStart	<- ifelse(bsubRelStart > 0.5, 1 - bsubRelStart, bsubRelStart)

# 	return(geneData)
# })
# sigF_allInAG_plot_df		<- bind_rows(sigF_allInAG_plot_list)



# Compared across replichores / All present / Some present / All absent genes
sigF_regulon_crossAG_comparison_plot	<- ggplot(subset(sigF_inAG_plot_df), aes(x = xmin, xend = xmax, y = relGeneStart, yend = relGeneStart)) +
	scale_x_continuous(name = "SigmaF regulon gene", limits = c(0, numGenes), breaks = seq(0.5, numGenes, by = 1), labels = unique(sigF_inAG_plot_df$geneName), sec.axis = dup_axis(name = NULL)) +
	scale_y_continuous(name = "Relative Genomic Position", limits = c(-0.05, 1.05), breaks = seq(0, 1, by = 0.5), labels = c("Origin", "Terminus", "Origin"), expand = c(0, 0)) +
	# Boundary lines
	geom_hline(yintercept = zoneBoundary_noPad_df$boundary, color = alpha(boundaryCol, 0.5), size = 0.5, linetype = "dashed") +
	# Zone colouring
	geom_rect(data = zoneBoundary_noPad_df, aes(xmin = -Inf, xmax = Inf, ymin = zoneMin, ymax = zoneMax), fill = alpha(zone_cols[1:9], 0.15), inherit.aes = FALSE) +
	# Gene positions
	geom_segment(color = alpha(wes_palette("Moonrise1")[4], 0.6)) +
	# Show which replichore the gene is on for Bsub
	geom_segment(
		data = subset(bSubOnly_pos_df, bsubReplichore == "Right"),
		aes(x = xmin + 0.4, xend = xmin + 0.4, y = -0.05, yend = 0.5, color = bsubReplichore),
		size = 0.3,
		inherit.aes = FALSE) +
	geom_segment(
		data = subset(bSubOnly_pos_df, bsubReplichore == "Left"),
		aes(x = xmin + 0.4, xend = xmin + 0.4, y = 0.5, yend = 1.05, color = bsubReplichore),
		size = 0.3,
		inherit.aes = FALSE) +
	coord_flip() +
	scale_color_manual(values = wes_palette("Darjeeling1")[c(1, 5)]) +
	theme_classic() +
	theme(
		panel.grid.major.y = element_line(size = 0.05),
		axis.text = element_text(size = 11),
		axis.title = element_text(size = 14))


# G. kaustophlius vs B. subtilis

gKauOnly_data	<- subset(sigF_inAG_plot_df, binomial == binomial_list[1])
gKauOnly_data$relStart	<- gKauOnly_data$relGeneStart

bSubOnly_data	<- subset(bSubOnly_pos_df, select = c(geneName, inAllAG, bsubRelStart, bsubToOriStart, bsubReplichore))
bSubOnly_data$binomial	<- "Bacillus subtilis 168"
bSubOnly_data$relStart	<- bSubOnly_data$bsubRelStart

gKay_bSub_comb		<- full_join(gKauOnly_data, subset(bSubOnly_data, select = c(binomial, geneName, inAllAG, bsubRelStart, bsubToOriStart, bsubReplichore, relStart)))

sigF_bSub_gKau_comparison_plot	<- ggplot(gKay_bSub_comb, aes(x = binomial, y = relStart, yend = relStart, group = geneName, color = binomial)) +
	scale_y_continuous(name = "Relative Genomic Position", limits = c(0, 1), breaks = seq(0, 1, by = 0.5), labels = c("Origin", "Terminus", "Origin")) +
	geom_point(shape = 124, size = 10) +
	geom_line(linetype = "dashed", color = alpha(wes_palette("IsleofDogs1")[6], 0.5)) +
	geom_label_repel(data = subset(gKay_bSub_comb, binomial == "Bacillus subtilis 168"),
		aes(x = 0.95, label = geneName),
		xlim = c(0, 0.7),
		hjust = 0,
		box.padding = 0.3,
		size = 3,
		show.legend = FALSE) +
	scale_color_manual(values = wes_palette("Darjeeling1")[c(1,5)]) +
	coord_flip() +
	theme_classic() +
	theme(
		panel.grid.major.y = element_line(size = 0.05),
		axis.text = element_text(size = 11),
		axis.title = element_text(size = 14))












ggplot(subset(sigF_inAG_plot_df), aes(x = xmin, xend = xmax, y = toOriStart, yend = toOriStart)) +
	scale_x_continuous(name = "SigmaF regulon gene", limits = c(0, numGenes), breaks = seq(0.5, numGenes, by = 1), labels = unique(sigF_inAG_plot_df$geneName)) +
	scale_y_continuous(name = "Relative Genomic Position", limits = c(0, 0.5), breaks = seq(0, 0.5, by = 0.5), labels = c("Origin", "Terminus")) +
	# Boundary lines
	geom_hline(yintercept = zoneBoundary_toOri_df$boundary, color = boundaryCol, size) +
	# Boundary padding
	geom_rect(data = zoneBoundary_toOri_df, aes(xmin = -Inf, xmax = Inf, ymin = boundaryMin, ymax = boundaryMax), fill = alpha(boundaryCol, 0.2), inherit.aes = FALSE) +
	# Zone colouring
	geom_rect(data = zoneBoundary_toOri_df, aes(xmin = -Inf, xmax = Inf, ymin = zoneMin, ymax = zoneMax), fill = zoneToOri_cols, inherit.aes = FALSE) +
	# Gene positions
	geom_segment(color = alpha(wes_palette("Moonrise1")[4], 0.6)) +
	coord_flip() +
	theme_classic() +
	theme(
		panel.grid.major.y = element_line(size = 0.05),
		axis.text = element_text(size = 11),
		axis.title = element_text(size = 14))













## Plot b subtilis genes only

bSubOnly_pos_list	<- lapply(unique(sigF_inAG_plot_df$geneName), function(gene) {
	perGene_data	<- subset(sigF_inAG_plot_df, geneName == gene)
	singleEntry		<- perGene_data[1,]
	return(singleEntry)
})
bSubOnly_pos_df		<- bind_rows(bSubOnly_pos_list)

ggplot(bSubOnly_pos_df, aes(x = 1, xend = 2, y = bsubToOriStart, yend = bsubToOriStart)) +
	scale_x_continuous(limits = c(1, 2.5)) +
	scale_y_reverse(limits = c(0.5, 0), breaks = seq(0.5, 0, by = -0.5), labels = c("Terminus", "Origin")) +
	geom_segment(color = alpha(wes_palette("Zissou1")[5], 0.6)) +
	geom_text_repel(
		aes(x = 2, y = bsubToOriStart, label = geneName),
		hjust = 0,
		direction = "y",
		nudge_x = 0.15,
		segment.size = 0.2,
		size = 2.5,
		box.padding = 0.6,
		max.iter = 10000) +
	theme_classic()


























