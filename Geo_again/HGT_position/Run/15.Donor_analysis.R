#!/usr/bin/env Rscript

# Load master variables and HGT position functions
invisible(sapply(HGTPos.all, source, .GlobalEnv))

# Any libraries explicitly used in the script
require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load("RSQLite", "tidyverse", "ggtree", "wesanderson", "Biostrings", "ggpubr", "Hmisc", "data.table", "broom")


# ------------------------------------------------------------------------------------- #
# Read in data

message("\nReading in data...", appendLF = FALSE)

# General position data
perTypeData			<- readRDS(file.path(positionData_path, "AG_perTypeData.rds"))

# Nodes and names NCBI taxon data
nodes_dt	<- readRDS(file.path(taxdmp_path, "nodes.rds"))
names_dt	<- readRDS(file.path(taxdmp_path, "names.rds"))

# Mowgli species tree
mowSpeciesTree_file	<- file.path(mowgli_path, "Mowgli_output", "Output_3", "1", "outputSpeciesTree.mpr")
if (!file.exists(mowSpeciesTree_file)) stop("Cannot find Mowgli species tree file, try a different file")
mowSpecies_tree		<- read.tree(mowSpeciesTree_file)
mowSpecies_phylo4	<- phylo4(mowSpecies_tree)

message("\rReading in data... done\n")

# ------------------------------------------------------------------------------------- #

uniqueEvents	<- unique(perTypeData$lHGT$'4'$allPosData$eventIndex)
donorBranches	<- perTypeData$lHGT$'4'$allPosData$donorEdge

donorID_list	<- mclapply(uniqueEvents, function(HGT_event) {


	donorBranch	<- perTypeData$lHGT$'4'$allPosData %>% filter(eventIndex == HGT_event) %>% distinct(donorEdge)

	branchNodes	<- unlist(str_split(donorBranch, " "))
	donorNode	<- branchNodes[2]

	allChildren <- descendants(mowSpecies_phylo4, donorNode)
	if (length(allChildren) == 0) {
		allChildren	<- str_subset(mowSpecies_tree$tip.label, paste0("_", donorNode))
		childTaxids	<- str_extract(allChildren, "[0-9]+")
	} else {
		childTaxids	<- str_extract(names(allChildren), "[0-9]+")
	}
	

	famNames_l	<- lapply(childTaxids, function(taxid) {
		taxidNum <- as.numeric(taxid)

		while (TRUE) {
			if (identical(taxidNum, 1)) break

			entry	<- nodes_dt %>% filter(taxid == taxidNum)
			rank	<- entry$nodeRank
			if (identical(rank, "family")) break
			taxidNum	<- entry$parentTaxid
		}

		if (identical(taxidNum, 1)) return(NA)

		familyTaxid	<- entry$taxid
		familyName	<- names_dt %>% filter(taxid == familyTaxid & nameClass == "scientific name") %>% pull(name)

		return(familyName)
	})

	famNames_l	<- famNames_l[!is.na(famNames_l)]
	if (length(famNames_l) == 0) {
		topFam <- NA
	} else {
		famName_tbl	<- sort(table(unlist(famNames_l)), decreasing = T)
		topFam		<- names(famName_tbl)[1]
		
	}
	out_df		<- data.frame(HGT_event = HGT_event, Family = topFam, stringsAsFactors = FALSE)
	return(out_df)
}, mc.cores = 20)


donorID_df	<- bind_rows(donorID_list)

commonDonors_tbl	<- sort(table(donorID_df$Family), decreasing = TRUE)[1:10]
commonDonors_df		<- as.data.frame(commonDonors_tbl)
names(commonDonors_df)	<- c("Family", "HGT_events")

commonDonorsBacillac_df	<- rbind(commonDonors_df, data.frame(Family = "Bacillaceae", HGT_events = 1608))
commonDonorsBacillac_df$Family	<- factor(commonDonorsBacillac_df$Family, levels = c("Bacillaceae", levels(commonDonors_df$Family)))

donorPlot <- ggplot(data = commonDonorsBacillac_df, aes(x = Family, y = HGT_events)) + geom_bar(stat = "identity") + lightTheme

quartz(width = 12, height = 8)
print(donorPlot)
quartz.save(file = "/Users/aesin/Desktop/Thesis/CH3/Figs/Figure_3H.pdf", type = "pdf", dpi = 300)
invisible(dev.off())


















