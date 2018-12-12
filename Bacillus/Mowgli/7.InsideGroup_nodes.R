#!/usr/local/bin/Rscript

require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load(tidyverse, fs, ape, phylobase)

## Main working dir
master			<- file.path("/Users/aesin/Desktop/Bacillus")

## Get the inside group AccAss list
insideGrp_dir		<- file.path(master, "Prefilter", "Bacillaceae")
insideAccAss_file	<- file.path(insideGrp_dir, "Bacillaceae_acc_ass_list.txt")
insideAccAss_tbl	<- read_tsv(file = insideAccAss_file, col_names = c("Acc_ass"))

## Get the Anoxy / Geobacillus taxid list
genomeLists_dir		<- file.path(master, "Core_genomes", "Genome_lists")
coreTaxid_file		<- file.path(genomeLists_dir, "coreToKeep.tsv")
coreTaxid_tbl		<- read_tsv(file = coreTaxid_file)

## Get the Acc_Ass & taxid translation table
accAssTransl_file	<- file.path(genomeLists_dir, "core_AccAssTaxid_table.tsv")
accAssTransl_tbl	<- read_tsv(accAssTransl_file)

## Using translation table, get the Inside Group taxids
insideTaxid_tbl		<- insideAccAss_tbl %>% left_join(accAssTransl_tbl, by = "Acc_ass")
## For Bacillus - remove any inside group nodes corresponding to trimmed Bacillus species (from within core group)
insideTaxid_tbl		<- insideTaxid_tbl %>% filter(!is.na(Taxid))

## Get the species tree Mowgli uses (any will do)
mowgliSpecies_file	<- file.path(master, "Mowgli", "Mowgli_output", "Output_3", "1", "outputSpeciesTree.mpr")
if (file.exists(mowgliSpecies_file)) {
	mowgliSpecies_tree	<- read.tree(mowgliSpecies_file)
} else {
	stop("Cannot find Mowgli species tree file, try a different file")
}

## Get the corresponding mowgli species tree tips (slightly different labelling format)
mowgliTips_tbl		<- tibble(mowTips = mowgliSpecies_tree$tip.label) %>%
	mutate(Taxid = str_extract(mowTips, "^[0-9]+(?=_)")) %>%
	mutate(Taxid = as.integer(Taxid))

insideTaxid_tbl		<- insideTaxid_tbl %>% left_join(mowgliTips_tbl, by = "Taxid")
coreTaxid_tbl		<- coreTaxid_tbl %>% left_join(mowgliTips_tbl, by = "Taxid")

# IGRenamedTips_list	<- unlist(lapply(insideTaxid_list, function(taxid) grep(pattern = paste0("^", taxid, "_"), mowgliTips_list, value = TRUE)))
# AGRenamedTips_list	<- unlist(lapply(AGtaxid_list, function(taxid) grep(pattern = paste0("^", taxid, "_"), mowgliTips_list, value = TRUE)))

## Prune just the inside group subtree and get all the node labels
insidePrune_tree	<- drop.tip(mowgliSpecies_tree, setdiff(mowgliSpecies_tree$tip.label, insideTaxid_tbl$mowTips))
insidePrune_as4		<- phylo4(insidePrune_tree)
insideNodeLab_tbl	<- tibble(insideNodes = as.numeric(nodeLabels(insidePrune_as4)))

## Make output directory
output_dir			<- file.path(direct, "Mowgli", "Inside_group")
dir_create(output_dir)

## Write inside group nodes out
write_delim(insideNodeLab_tbl, path = file.path(output_dir, "InsideGroup_nodes.txt"), col_names = FALSE, delim = "\n")
write_tsv(insideTaxid_tbl, path = file.path(output_dir, "InsideGroup_full.tsv"))
write_tsv(coreTaxid_tbl, path = file.path(output_dir, "Core_full.tsv"))