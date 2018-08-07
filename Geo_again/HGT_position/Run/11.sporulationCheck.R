#!/usr/bin/env Rscript

# Load master variables and HGT position functions
invisible(sapply(HGTPos.all, source, .GlobalEnv))

# Any libraries explicitly used in the script
require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load("dplyr", "reshape2", "grid", "ggplot2", "ggrepel", "ggdendro", "ggpubr", "wesanderson", "RColorBrewer")

# ------------------------------------------------------------------------------------- #
message("\nReading in data...", appendLF = FALSE)

# General position data
perTypeData			<- readRDS(file.path(positionData_path, "AG_perTypeData.rds"))

# Subgroup data
dnaA_pos_data		<- readRDS(file.path(positionData_path, "bySpecies_dnaA_data.rds"))

# Zone boundaries
zoneBoundaryList	<- readRDS(file.path(positionData_path, "AG_zoneBoundaries.rds"))

message("\rReading in data... done\n")
# ------------------------------------------------------------------------------------- #
# Open All_prot database
dbConn			<- dbConnect(RSQLite::SQLite(), allProtDB_path)

binomial_list		<- unique(perTypeData$All$allPosData$binomial)

# Path to sporulation analysis
sporFigureOut_path	<- file.path(sporulation_path, "Figures", "23Genomes")
if (!dir.exists(sporFigureOut_path)) dir.create(sporFigureOut_path, recursive = TRUE)


# ------------------------------------------------------------------------------------- #

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
	if (nrow(AG_inOrth) == 0) {
		emptyRow			<- subset(perTypeData$All$allPosData[1,], select = -c(CircStart, CircEnd))
		emptyRow[1,]		<- NA
		emptyRow$geneName	<- sigF_geneName
		emptyRow$inAllAG	<- FALSE
		return(list(perOrth_data = emptyRow, uniqueTaxa = 0))
	}

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

# ------------------------------------------------------------------------------------- #

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

# Compared across replichores / All present / Some present / All absent genes
sigF_regulon_crossAG_comparison_plot	<- ggplot(subset(sigF_inAG_plot_df), aes(x = xmin, xend = xmax, y = relGeneStart, yend = relGeneStart)) +
	scale_x_continuous(name = "SigmaF regulon gene", limits = c(0, numGenes), breaks = seq(0.5, numGenes, by = 1), labels = unique(sigF_inAG_plot_df$geneName), sec.axis = dup_axis(name = NULL)) +
	scale_y_continuous(name = "Relative Genomic Position", limits = c(-0.05, 1.05), breaks = seq(0, 1, by = 0.5), labels = c("Origin", "Terminus", "Origin"), expand = c(0, 0)) +
	# Boundary lines
	geom_hline(yintercept = zoneBoundaryList$fullRange$boundary, color = alpha(zoneBoundaryList$boundCol, 0.5), size = 0.5, linetype = "dashed") +
	# Zone colouring
	geom_rect(data = zoneBoundaryList$fullRange, aes(xmin = -Inf, xmax = Inf, ymin = zoneMin, ymax = zoneMax), fill = zoneBoundaryList$fullRange$zoneCol_alpha, inherit.aes = FALSE) +
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


# ------------------------------------------------------------------------------------- #
# G. kaustophlius vs B. subtilis

gKauOnly_sigF_data	<- subset(sigF_inAG_plot_df, binomial == binomial_list[1])
gKauOnly_sigF_data$relStart	<- gKauOnly_sigF_data$relGeneStart

bSubOnly_sigF_data	<- subset(bSubOnly_pos_df, select = c(geneName, inAllAG, bsubRelStart, bsubToOriStart, bsubReplichore))
bSubOnly_sigF_data$binomial	<- "Bacillus subtilis 168"
bSubOnly_sigF_data$relStart	<- bSubOnly_sigF_data$bsubRelStart

gKau_bSub_sigF_comb		<- full_join(gKauOnly_sigF_data, subset(bSubOnly_sigF_data, select = c(binomial, geneName, inAllAG, bsubRelStart, bsubToOriStart, bsubReplichore, relStart)))

sigF_bSub_gKau_comparison_plot	<- ggplot(gKau_bSub_sigF_comb, aes(x = binomial, y = relStart, yend = relStart, group = geneName, color = binomial)) +
	scale_y_continuous(name = "Relative Genomic Position", limits = c(0, 1), breaks = seq(0, 1, by = 0.5), labels = c("Origin", "Terminus", "Origin")) +
	geom_point(shape = 124, size = 10) +
	geom_line(linetype = "dashed", color = alpha(wes_palette("IsleofDogs1")[6], 0.5)) +
	geom_label_repel(data = subset(gKau_bSub_sigF_comb, binomial == "Bacillus subtilis 168"),
		aes(x = 0.95, label = geneName),
		xlim = c(0, 0.7),
		hjust = 0,
		box.padding = 0.3,
		size = 3,
		show.legend = FALSE) +
	scale_color_manual(values = wes_palette("Darjeeling1")[c(1,5)]) +
	coord_flip() +
	ggtitle("Comparison of sigmaF regulon genes between Geobacillus kaustophilus and Bacillus subtilis") +
	theme_classic() +
	theme(
		panel.grid.major.y = element_line(size = 0.05),
		axis.text = element_text(size = 11),
		axis.title = element_text(size = 14),
		plot.title = element_text(hjust = 0.5, size = 16))



# ------------------------------------------------------------------------------------- #
# There are 1098 vertical genes in Gkau - corresponding to 1094 orthGroups
gKauOnly_ver_data	<- subset(perTypeData$Ver$'3'$allPosData, binomial == binomial_list[1])
gKauOnly_verGroup_l	<- unique(gKauOnly_ver_data$orthGroup)

# A full list of ortholog genes to write out
bsub_allVer_protID_tbl	<- dbSendQuery(dbConn, 'SELECT * FROM t1 WHERE taxid == :bsub_taxid AND OrthGroup == :orthGroup')
dbBind(bsub_allVer_protID_tbl, param = list(bsub_taxid = rep(bsub_168_taxid, length(gKauOnly_verGroup_l)), orthGroup = gKauOnly_verGroup_l))
bsub_allVer_raw_df	<- dbFetch(bsub_allVer_protID_tbl) %>% select(-c(sequence, NuclSeq))
dbClearResult(bsub_allVer_protID_tbl)

# ---------------------------------- #

# Find the relative start position of the B subtilis genes
bsub_ver_relStart_list	<- lapply(bsub_allVer_raw_df$gene_start, genomeRelativePosition, oriStart = bsub_168_dnaA_data$oriStart, oriEnd = bsub_168_dnaA_data$oriEnd, oriStrand = bsub_168_dnaA_data$oriStrand, genomeLength = bsub_168_genomeLength)

# Adjust the dataframe to conform to the gKau data
bsub_allVer_adj_df	<- bsub_allVer_raw_df %>%
	mutate(relStart = unlist(bsub_ver_relStart_list)) %>%
	mutate(distToOri = case_when(
		relStart > 0.5 ~ 1 - relStart,
		relStart <= 0.5 ~ relStart)) %>%
	mutate(binomial = "Bacillus subtilis 168") %>%
	select(-c(is_ag, strain, acc_ass, product)) %>%
	dplyr::rename(orthGroup = OrthGroup, geneStart = gene_start, geneEnd = gene_end, locusTag = locus)


# Save as RDS
write.csv(bsub_allVer_adj_df, file = file.path(sporulation_path, "bSubtilus_fullOrth_list.csv"), row.names = FALSE)

# ---------------------------------- #
# For 1094 Gkau Vertical orthGroups there are 898 Bsub genes from 893 groups. Select groups only present in both
gKauOnly_verInBsub_data				<- subset(gKauOnly_ver_data, orthGroup %in% bsub_allVer_adj_df$orthGroup)
gKauOnly_verInBsub_data$relStart	<- gKauOnly_verInBsub_data$relGeneStart

# Join the Gkau and Bsub data
gKau_bSub_allVer_comb		<- full_join(gKauOnly_verInBsub_data, bsub_allVer_adj_df %>% select(-COGcat))

# ---------------------------------- #
# Find gKau genome length
gKau_genomeLength_tbl	<- dbSendQuery(dbConn, 'SELECT genome_l FROM t1 WHERE protId == :sampleProtID LIMIT 1')
dbBind(gKau_genomeLength_tbl, param = list(sampleProtID = sample(gKauOnly_ver_data$protID, 1)))
gKau_genomeLength		<- dbFetch(gKau_genomeLength_tbl) %>% pull(genome_l)
dbClearResult(gKau_genomeLength_tbl)


# Scale by distance from terminus
gKau_bSub_final	<- gKau_bSub_allVer_comb %>% 
	mutate(genome_l = replace(genome_l, is.na(genome_l), gKau_genomeLength)) %>%
	mutate(terRelStart = geneStart - (genome_l / 2))

bsub_sigF_relTer	<- bsub_sigF_full_df %>%
	mutate(terRelStart = gene_start - (bsub_168_genomeLength / 2)) %>%
	mutate(binomial = "Bacillus subtilis 168") %>%
	dplyr::rename(orthGroup = OrthGroup)

# Plot
quartz(width = 21, height = 8, canvas = "white", bg = "white")
allVer_bSub_gKau_comparison_plot	<- ggplot(gKau_bSub_final, aes(x = binomial, y = terRelStart, group = orthGroup, color = binomial)) +
	scale_y_continuous(
		name = "Position Relative to Terminus") +
	geom_point(shape = 124, size = 8) +
	geom_point(data = bsub_sigF_relTer, aes(x = 0.5, y = terRelStart, group = orthGroup), shape = 124, size = 8, color = "orange") +
	geom_line(linetype = "solid", color = alpha(wes_palette("IsleofDogs1")[6], 0.9), size = 0.1) +
	scale_color_manual(values = wes_palette("Darjeeling1")[c(5,1)], guide = FALSE) +
	coord_flip() +
	ggtitle("Comparison of shared (AG vertical) genes between Geobacillus kaustophilus and Bacillus subtilis") +
	theme_classic() +
	theme(
		panel.grid.major.y = element_line(size = 0.05),
		panel.grid.major.x  = element_line(size = 0.4, linetype = "dashed", color = wes_palette("Darjeeling1")[2]),
		axis.text = element_text(size = 11),
		axis.text.y = element_text(angle = 45),
		axis.title = element_text(size = 14),
		axis.title.y = element_blank(),
		plot.title = element_text(hjust = 0.5, size = 16))
print(allVer_bSub_gKau_comparison_plot)

fileName	<- file.path(sporFigureOut_path, "allVer_bSub_gKau_terCentred")
invisible(quartz.save(file = paste0(fileName, ".pdf"), type = "pdf", dpi = 100))
invisible(quartz.save(file = paste0(fileName, ".png"), type = "png", dpi = 100))
invisible(dev.off())


# ------------------------------------------------------------------------------------- #
# Bacillus subtilus expression analysis
geneExpr_data	<- read.delim(file.path(sporulation_path, "GSE78108_growth_rate_annotated_genes.txt"), header = TRUE, stringsAsFactors = FALSE)
allGene_data	<- read.delim(file.path(sporulation_path, "All_genes_of_B._subtilis_subtilis_168.txt"), header = TRUE, stringsAsFactors = FALSE)

# Rename to locus tag
allGene_data	%<>% dplyr::rename(Locus.Tag = Accession.1) %>% arrange(Left.End.Position)

# Combine expression data with bsub position data
combined_data		<- left_join(geneExpr_data, allGene_data[,2:4], by = "Locus.Tag")
combined_data$CHG	<- log(rowSums(combined_data[,grep("CHG\\.medium", names(combined_data))]) / 4)
combined_data$CH	<- log(rowSums(combined_data[,grep("CH\\.medium", names(combined_data))]) / 4)
combined_data$S		<- log(rowSums(combined_data[,grep("S\\.medium", names(combined_data))]) / 4)
combined_data$SE	<- log(rowSums(combined_data[,grep("SE\\.medium", names(combined_data))]) / 4)

# Get genome length and normalise genome size
combined_data		%<>% mutate(relGeneStart = Left.End.Position / bsub_168_genomeLength)

# Extend the data to avoid edge effects
combined_extend		<- combined_data %>%
	subset(relGeneStart <= 0.125 | relGeneStart >= 0.875) %>%
	mutate(relGeneStart = case_when(
		relGeneStart < 0.5 ~ relGeneStart + 1,
		relGeneStart > 0.5 ~ relGeneStart - 1)) %>%
	bind_rows(combined_data, .) %>%
	arrange(relGeneStart)

# Take our Bsub vertical genes and find the corresponding gene in the expression dataset by checking the closest genes
allGene_data_starts	<- allGene_data$Left.End.Position
num_allGeneData		<- nrow(allGene_data)

withAdjustedGenePositions	<- lapply(1:nrow(bsub_allVer_raw_df), function(geneIndex) {
	geneEntry	<- bsub_allVer_raw_df[geneIndex,]
	geneStart	<- geneEntry$gene_start

	absDists	<- abs(rep(geneStart, num_allGeneData) - allGene_data_starts)
	closestGene	<- allGene_data[which(absDists == min(absDists)),]

	geneEntry$adj_distance	<- closestGene$Left.End.Position
	geneEntry$relAdj_dist	<- genomeRelativePosition(genePosition = geneEntry$adj_distance, oriStart = bsub_168_dnaA_data$oriStart, oriEnd = bsub_168_dnaA_data$oriEnd, oriStrand = bsub_168_dnaA_data$oriStrand, genomeLength = bsub_168_genomeLength)

	return(geneEntry)
})
withAdjustedGenePositions_df	<- bind_rows(withAdjustedGenePositions)

combined_extend.melt	<- melt(combined_extend, id.vars = c("Gene.Name", "relGeneStart"), measure.vars = c("CHG", "CH", "S", "SE"))

plotCols	<- brewer.pal(6, "Blues")
bSub_expression_plot	<- ggplot(data = combined_extend.melt, aes(x = relGeneStart, y = value, color = variable)) +
	scale_y_continuous(name = "Relative Gene Expression") +
	geom_smooth(method = "loess", span = 0.2, se = FALSE) +
	geom_point(data = withAdjustedGenePositions_df, aes(x = relAdj_dist, y = 7), shape = 124, size = 8, inherit.aes = FALSE) +
	scale_color_manual(name = "Growth\nMedium", values = rev(plotCols[3:6])) +
	coord_cartesian(xlim = c(0, 1), ylim = c(4.9, 9), expand = FALSE) +
	lightTheme +
	theme(
		panel.grid.minor.y = element_blank(),
		panel.grid.major.x  = element_line(),
		axis.text = element_text(size = 11),
		axis.title = element_text(size = 14),
		axis.title.y = element_blank(),
		plot.title = element_text(hjust = 0.5, size = 16)
	)

quartz(width = 21, height = 8, canvas = "white", bg = "white")
print(bSub_expression_plot)

fileName	<- file.path(sporFigureOut_path, "bSub_expression_plot")
invisible(quartz.save(file = paste0(fileName, ".pdf"), type = "pdf", dpi = 100))
invisible(quartz.save(file = paste0(fileName, ".png"), type = "png", dpi = 100))
invisible(dev.off())













