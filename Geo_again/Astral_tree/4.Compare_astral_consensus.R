#!/usr/local/bin/Rscript

library(ape)
library(phytools)
library(phylobase)
library(phangorn)
library(stringr)

NumericTaxidString	<- function(taxid_vector) {
	new_vector	<- paste("tax", taxid_vector, sep = "_")
	return(new_vector)
}

StringTaxidNumeric	<- function(taxid_ST_vector) {
	new_vector	<- str_sub(taxid_ST_vector, 5, -1)
	return(new_vector)
}

direct			<- "/Users/aesin/Desktop/Geo_again"
spec_tree_taxs	<- file.path(direct, "Astral_speciesTree", "Astral_noDup_speciesTree_ultra.txt")
bacill_tree_dir	<- file.path(direct, "Consensus_groups", "Bacillaceae", "Final_trees", "Taxid_labelled")
output_dir		<- file.path(direct, "Astral_speciesTree", "Correct_anogeo")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

## Get a list of the Bacillaceae taxids
# Bacillaceae acc_ass and the translation table files
bacillac_acc_ass	<- file.path(direct, "Consensus_groups", "Bacillaceae", "Bacillaceae_acc_ass_list.txt")
anogeo_taxid_file	<- file.path(direct, "Genomes", "Genome_lists", "AG_taxids.txt")
acc_ass_taxid_trans	<- file.path(direct, "Genomes", "Genome_lists", "Acc_ass_taxid_table.tsv")
taxid_binom_trans	<- file.path(direct, "Genomes", "Genome_lists", "Taxid_refinedBinomial_table.tsv")

# Read in the acc_ass and translation files
bac_acc_ass_list	<- scan(bacillac_acc_ass, what = character())
acc_ass_transl		<- read.table(acc_ass_taxid_trans, sep = "\t", header = FALSE, col.names = c("Acc_ass", "Taxid"), colClasses = rep("character", 2))
taxid_bin_transl	<- read.table(taxid_binom_trans, sep = "\t", header = TRUE, colClasses = rep("character", 2))

# Use the acc_ass to get the Bacillaceae taxids
bacillac_taxids		<- acc_ass_transl$Taxid[which(acc_ass_transl$Acc_ass %in% bac_acc_ass_list)]
bacillac_taxids_ST	<- NumericTaxidString(bacillac_taxids)

## Read in Astral species tree and get the subtree containing all Bacillaceae
astral_taxs_data			<- read.tree(spec_tree_taxs)
astral_taxs_data$tip.label	<- NumericTaxidString(astral_taxs_data$tip.label)
# Use the phytools midpoint.root - the phangorn midpoint() produces a weird tree
astral_taxs_MP				<- midpoint.root(astral_taxs_data)

## Subset all the taxa from the Astral tree corresponding to the Bacillaceae
## use the most recent common ancestor to see how compactly the Bacillaceae
## are represented in the Astral data.
astral_taxs_MP4		<- phylo4(astral_taxs_MP)
bacillacMRCA		<- MRCA(astral_taxs_MP4, bacillac_taxids_ST)
subset_tree 		<- subset(astral_taxs_MP4, mrca = bacillacMRCA)
bacillac_only		<- as(subset(subset_tree, tips.include = bacillac_taxids_ST, trim.internal = TRUE), "phylo")

## The Bacillaceae are highly associated within the Astral tree - the 215 tips
## are represented fully in a subset tree of 236 tips - i.e. the Bacillaceae are
## not quite monophyletic but the extra taxa pulled out with them appears to also
## be closely related
subset_tips			<- StringTaxidNumeric(as(subset_tree, "phylo")$tip.label)
extra_tips			<- setdiff(subset_tips, bacillac_taxids)
extra_tips_binoms	<- taxid_bin_transl$Binomial[which(taxid_bin_transl$Taxid %in% extra_tips)]

## Let's check the similarity of the monophyletic AnoGeo branching in the Astral
## tree to the branching in the RAxML-reconstructed Bacillaceae tree(s)
anogeo_taxids		<- scan(anogeo_taxid_file, what = character())
anogeo_taxids_ST	<- NumericTaxidString(anogeo_taxids)

tree_evalues		<- list("-10", "-50", "-100", "-150")
consensus_trees	<- lapply(tree_evalues, function(trunc_eval) {

	tree_type		<- list("RAxML", "FT")
	per_type		<- lapply(tree_type, function(type) {
		if (type == "RAxML") {
			tree_file		<- file.path(bacill_tree_dir, paste0(type, "_bipartitions.super_tree", trunc_eval, ".txt"))	
		} else {
			tree_file		<- file.path(bacill_tree_dir, paste0(type, "_super_tree", trunc_eval, ".txt"))	
		}
		tree_data		<- read.tree(tree_file)
		tree_data_MP	<- midpoint.root(tree_data)
		if (is.binary.phylo(tree_data_MP) == FALSE) {
			tree_data_MP	<- multi2di(tree_data_MP)
		}
		## Relabel the tips so that phylo4 doesnt get confused
		tree_data_MP$tip.label	<- NumericTaxidString(tree_data_MP$tip.label)
		## Drop the bootstrap values as node labels
		tree_data_4		<- phylo4(tree_data_MP, check.node.labels = "drop")
		anogeoMRCA		<- MRCA(tree_data_4, anogeo_taxids_ST)
		anogeo_sub		<- subset(tree_data_4, mrca = anogeoMRCA)
		anogeo_phylo	<- as(anogeo_sub, "phylo")
		return(list(bacillac = tree_data_MP, anogeo = anogeo_phylo))
	})
	names(per_type) <- tree_type
	return(per_type)
})
names(consensus_trees)	<- tree_evalues

## RF distances between bacillaceae consensus trees for the same evalue
RF_anogeo_raxft_50	<- RF.dist(consensus_trees$`-10`$FT$anogeo, consensus_trees$`-10`$RAxML$anogeo)
RF_anogeo_raxft_50	<- RF.dist(consensus_trees$`-50`$FT$anogeo, consensus_trees$`-50`$RAxML$anogeo)
RF_anogeo_raxft_100	<- RF.dist(consensus_trees$`-100`$FT$anogeo, consensus_trees$`-100`$RAxML$anogeo)
RF_anogeo_raxft_150	<- RF.dist(consensus_trees$`-150`$FT$anogeo, consensus_trees$`-150`$RAxML$anogeo)

## Extract anogeo from the Astral tree
anogeoMRCA			<- MRCA(astral_taxs_MP4, anogeo_taxids_ST)
anogeo_subset_tree 	<- subset(astral_taxs_MP4, mrca = anogeoMRCA)
anogeo_subset_phylo	<- as(anogeo_subset_tree, "phylo")

RAxML_anogeo_trees	<- list(consensus_trees$`-10`$RAxML$anogeo, consensus_trees$`-50`$RAxML$anogeo, consensus_trees$`-100`$RAxML$anogeo, consensus_trees$`-150`$RAxML$anogeo, anogeo_subset_phylo)
names(RAxML_anogeo_trees)	<- c(tree_evalues, "astral")
class(RAxML_anogeo_trees)	<- "multiPhylo"
cross_RF_dists		<- RF.dist(RAxML_anogeo_trees)

## We see that the only difference between -10 + -50 + Astral Anogeo phylogenies
## is the exact position of G. GHH01 - with no consensus, we got to the fuller
## anogeo phylogeny reconstructions (using non-assembled genomes as well). There
## we see a consensus of G. GHH01 branching above both the G.kay/G.ther and
## the G.C56-T3 clades - corresponding to the topolgy in the -50 consensus tree above.

## We correct the Astral tree Anogeo phylogeny to correspond to the -50 topology.
astral_phylo		<- as(astral_taxs_MP4, "phylo")	
prune_targets		<- anogeo_subset_phylo$tip.label[1:2]
almost_pruned		<- drop.tip(astral_phylo, anogeo_subset_phylo$tip.label[-1:-2])
almost_pruned$tip.label[which(almost_pruned$tip.label %in% prune_targets)] <- c("BIND1", "BIND2")

rebind_where		<- getMRCA(almost_pruned, c("BIND1", "BIND2"))
rebind_anogeo		<- bind.tree(almost_pruned, consensus_trees$`-50`$RAxML$anogeo, where = rebind_where, position = 0.1)
rebind_anogeo		<- multi2di(rebind_anogeo, random = FALSE)

rebind_prune_keys	<- drop.tip(rebind_anogeo, c("BIND1", "BIND2"))
rebind_prune_keys$tip.label	<- StringTaxidNumeric(rebind_prune_keys$tip.label)
# write.tree(rebind_prune_keys, file = file.path(output_dir, "Astral_speciesTree_corrected.txt"))

## After writing out - put the file through FigTree to make the tree ultrametric
## Below is optional code to read that tree in and check it's ultrametric

# ultra_correct_tree	<- read.tree(file = file.path(output_dir, "Astral_speciesTree_correct_ultra.txt"))
# is.ultrametric(ultra_correct_tree)

# Visually very few disagreements
x <- cophylo(bacillac_only, consensus_trees$`-50`$RAxML$bacillac, rotate = FALSE)

quartz(h = 40, w = 10)
plot.cophylo(x, size = 2)
quartz.save(file = file.path("/Users/aesin/Desktop/testCophylo.pdf"), type = "pdf", dpi = 300)

# Will need to thoroughly check whuch groups are included in the Astral bacillaceae MRCA and are not bacillaceae
# The 21 extra tips appear to all be Bacillales (need to double check with lit)













