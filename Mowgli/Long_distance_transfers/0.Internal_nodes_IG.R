library(ape)
library(phylobase)
library(stringr)

#### WRITE OUT THE INTERNAL NODE LABELS FOR THE INSIDE GROUP. If the file Inside_group_nodes.txt already exists, does not need to be re-run ####

setwd("/Users/aesin/Desktop/Mowgli/")

## If the output node file exists, abort the script ##
if (file.exists("Inside_group_nodes.txt") == TRUE) stop("Output file \"Inside_group_nodes.txt\" already exists. Exiting.")

## Read in the list of inside group tips and create an empty vector ##
inside_group_tips	<- as.vector(read.delim("Inside_group_tips.txt", header = FALSE)$V1)
renamed_IG_tips		<- vector()

## Read in the species tree. We have to use a species tree processed by reconciliation (this assigns the node labels) ##
species_tree <- read.tree(file = "Species_tree/Reconciled_species_tree/outputSpeciesTree.mpr")

## Tip names were changed for the species tree to be used by Mowgli. But, the original name is always a subset of the name used in the species tree, so we can find the list of modified labels with grep ##
for (tip in inside_group_tips) {

	renamed_tip		<- grep(tip, species_tree$tip.label, value = TRUE)
	renamed_IG_tips	<- c(renamed_IG_tips, renamed_tip)
}

## Using this list of tip labels, we can get a prune of all the IG group by dropping all non-IG tips and get all the node labels from that ##
IG_prune	<- drop.tip(species_tree, setdiff(species_tree$tip.label, renamed_IG_tips))
IG_prune_4	<- phylo4(IG_prune)
node_labels	<- as.vector(as.numeric(nodeLabels(IG_prune_4)))

## Write the node labels out ##
write.table(node_labels, file = "Inside_group_nodes.txt", quote = FALSE, row.names = FALSE, sep = "\n")
