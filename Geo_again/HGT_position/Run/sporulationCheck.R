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
bsub_sigF_full_df$relStart	<- unlist(bsub_sigF_relStart_list)

# Sample plot for G sub (compare to Wang JMB 2006 Fig 1)
ggplot(bsub_sigF_full_df, aes(x = relStart, xend = relStart, y = 1.2, yend = 1, label = bSub_sigF_genes)) + geom_segment() + coord_polar() + scale_y_continuous(limits = c(0,1.25)) + geom_text_repel(size = 3)

# ------------------------------------------------------------------------------------- #
# Of the 48 genes from Wang, 43 are represented in Geobacillus orthologous groups
# The 43 genes belong to seperate orthologous families in AG

sigF_gene_missingAG	<- bsub_sigF_full_df[is.na(bsub_sigF_full_df$OrthGroup),]
sigF_gene_presentAG	<- bsub_sigF_full_df[!is.na(bsub_sigF_full_df$OrthGroup),]

# In 32 of the 43 orthologous groups, all taxa were represented
sigF_inAG_data <- lapply(1:nrow(sigF_gene_presentAG), function(sigF_index) {
	sigF_gene	<- sigF_gene_presentAG[sigF_index,]

	sigF_orthGroup	<- sigF_gene$OrthGroup
	sigF_geneName	<- sigF_gene$bSub_sigF_genes

	AG_inOrth	<- subset(perTypeData$All$allPosData, orthGroup == sigF_orthGroup & !taxid %in% outlierTaxid, select = -c(CircStart, CircEnd))
	AG_inOrth$geneName	<- sigF_geneName

	uniqueTaxa	<- length(unique(AG_inOrth$binomial))
	print(uniqueTaxa)
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

ggplot(subset(sigF_inAG_df, inAllAG == TRUE), aes(x = relGeneStart, xend = relGeneStart, y = 1.2, yend = 1)) + geom_segment() + coord_polar() + scale_y_continuous(limits = c(0,1.25))
