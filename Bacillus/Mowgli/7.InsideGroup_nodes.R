#!/usr/local/bin/Rscript

if (!require("pacman")) install.packages("pacman")
pacman::p_load("ape", "stringr", "phylobase")

## Main working dir
direct			<- file.path("/Users/aesin/Desktop/Geo_again")

## Get the inside group AccAss list
insideGrp_dir		<- file.path(direct, "Consensus_groups", "Bacillaceae")
insideAccAss_file	<- file.path(insideGrp_dir, "Bacillaceae_acc_ass_list.txt")
insideAccAss_list	<- read.table(file = insideAccAss_file, sep = "\n", stringsAsFactors = FALSE)$V1

## Get the Anoxy / Geobacillus taxid list
genomeLists_dir		<- file.path(direct, "Genomes", "Genome_lists")
AGtaxid_file		<- file.path(genomeLists_dir, "AG_taxids.txt")
AGtaxid_list		<- read.table(file = AGtaxid_file, sep = "\n", stringsAsFactors = FALSE)$V1

## Get the Acc_Ass & taxid translation table
accAssTransl_file	<- file.path(genomeLists_dir, "Acc_ass_taxid_table.tsv")
accAssTransl_df		<- read.table(accAssTransl_file, sep = "\t", stringsAsFactors = FALSE, col.names = c("AccAss", "Taxid"))

## Using translation table, get the Inside Group taxids
insideTaxid_list	<- accAssTransl_df$Taxid[which(accAssTransl_df$AccAss %in% insideAccAss_list)]

## Get the species tree Mowgli uses (any will do)
mowgliSpecies_file	<- file.path(direct, "Mowgli", "Mowgli_output", "Output_3", "1", "outputSpeciesTree.mpr")
if (file.exists(mowgliSpecies_file)) {
	mowgliSpecies_tree	<- read.tree(mowgliSpecies_file)
} else {
	stop("Cannot find Mowgli species tree file, try a different file")
}

## Get the corresponding mowgli species tree tips (slightly different labelling format)
mowgliTips_list		<- mowgliSpecies_tree$tip.label
IGRenamedTips_list	<- unlist(lapply(insideTaxid_list, function(taxid) grep(pattern = paste0("^", taxid, "_"), mowgliTips_list, value = TRUE)))
AGRenamedTips_list	<- unlist(lapply(AGtaxid_list, function(taxid) grep(pattern = paste0("^", taxid, "_"), mowgliTips_list, value = TRUE)))

## Prune just the inside group subtree and get all the node labels
insidePrune_tree	<- drop.tip(mowgliSpecies_tree, setdiff(mowgliTips_list, IGRenamedTips_list))
insidePrune_as4		<- phylo4(insidePrune_tree)
insideNodeLab_list	<- as.vector(as.numeric(nodeLabels(insidePrune_as4)))

## Make output directory
output_dir			<- file.path(direct, "Mowgli", "Inside_group")
dir.create(output_dir, showWarnings = FALSE)

## Write inside group nodes out
write.table(insideNodeLab_list, file = file.path(output_dir, "InsideGroup_nodes.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\n")
write.table(IGRenamedTips_list, file = file.path(output_dir, "InsideGroup_mowTips.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\n")
write.table(AGRenamedTips_list, file = file.path(output_dir, "AnoxyGeo_mowTips.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\n")