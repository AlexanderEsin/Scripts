#!/usr/local/bin/Rscript

require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load(tidyverse, fs, ape, phylobase)

## Main working dir
master			<- file.path("/Users/aesin/Desktop/Staph")

## Get the inside group AccAss list
specTree_dir	<- file.path(master, "Astral_speciesTree")
specTree_file	<- file.path(specTree_dir, "Astral_speciesTree_ultraChecked.txt")
specTree_tree	<- read.tree(specTree_file)

## Get the Anoxy / Geobacillus taxid list
genomeLists_dir		<- file.path(master, "Core_genomes", "Genome_lists")
coreTaxid_file		<- file.path(genomeLists_dir, "coreToKeep.tsv")
coreTaxid_tbl		<- read_tsv(file = coreTaxid_file)

# Find the outgroups to the Staph species (macrococcus etc.) 1785995 is the most removed outgroup; 93061 is a staph species
staphWithOut_node	<- getMRCA(specTree_tree, c("1785995", "93061"))
staphWithOut_clade	<- extract.clade(specTree_tree, staphWithOut_node)
insideTaxids_tbl	<- tibble(Taxid = as.integer(staphWithOut_clade$tip.label))

## Get the Acc_Ass & taxid translation table
accAssTransl_file	<- file.path(genomeLists_dir, "core_AccAssTaxid_table.tsv")
accAssTransl_tbl	<- read_tsv(accAssTransl_file)

## Using translation table, get the Inside Group taxids
insideTaxid_tbl		<- insideTaxids_tbl %>% left_join(accAssTransl_tbl, by = "Taxid")

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

## Prune just the inside group subtree and get all the node labels
insidePrune_tree	<- drop.tip(mowgliSpecies_tree, setdiff(mowgliSpecies_tree$tip.label, insideTaxid_tbl$mowTips))
insidePrune_as4		<- phylo4(insidePrune_tree)
insideNodeLab_tbl	<- tibble(insideNodes = as.numeric(nodeLabels(insidePrune_as4)))

## Make output directory
output_dir			<- file.path(master, "Mowgli", "Inside_group")
dir_create(output_dir)

## Write inside group nodes out
write_delim(insideNodeLab_tbl, path = file.path(output_dir, "InsideGroup_nodes.txt"), col_names = FALSE, delim = "\n")
write_tsv(insideTaxid_tbl, path = file.path(output_dir, "InsideGroup_full.tsv"))