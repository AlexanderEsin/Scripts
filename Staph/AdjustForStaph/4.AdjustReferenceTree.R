#!/usr/bin/env Rscript

# Any libraries explicitly used in the script
require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load(tidyverse, fs, ape)

master			<- file.path("/Users/aesin/Desktop/Bacillus")
bac_list_dir	<- file.path(master, "Bac_genomes", "Genome_lists")

# Full (final) Astral species tree from the AG set
fullAstral_tree	<- file.path(master, "..", "Geo_again", "Astral_speciesTree", "Correct_anogeo", "Astral_speciesTree_correct_ultra.txt")
fullAstral_data	<- read.tree(fullAstral_tree)

# New Astral reference tree for this group reconciliation
newAstra_dir	<- file.path(master, "Astral_speciesTree")
dir_create(newAstra_dir)

# The bacillus consensus tree
finalBacil_tree	<- file.path(master, "Prefilter", "Bacillus_withOut", "Final_trees", "Taxid_labelled", "FT_super_tree-50.txt")
finalBacil_data	<- read.tree(finalBacil_tree)

# Tables to read in
toRemove			<- read_tsv(file.path(bac_list_dir, "bacillusToRemove.tsv"))
toKeep				<- read_tsv(file.path(bac_list_dir, "bacillusToKeep.tsv"))
bac_taxidBinom_tbl	<- read_tsv(file.path(bac_list_dir, "bac_TaxidBinomial_table.tsv"))

# ----------------------------------------------------- #
# First remove the excess tips that are no longer included in the dataset (extra bacillus)
astralPrune_data	<- drop.tip(fullAstral_data, tip = as.character(toRemove$Taxid))

# Check whether our bacillus is monophyletic & extract clude from Astral tree
message(paste0("Group is monophyletic in Astral: ", is.monophyletic(astralPrune_data, tip = as.character(toKeep$Taxid))))
astr_mrca	<- getMRCA(astralPrune_data, tip = as.character(toKeep$Taxid))
astr_clade	<- extract.clade(astralPrune_data, node = astr_mrca)

# Now extract the relevant tips from the Bacillus consensus tree
# First find the root node of our bacillus group (the consensus tree has na outgroup)
finalBacil_mp	<- midpoint(finalBacil_data)
bac_mrca		<- getMRCA(finalBacil_mp, tip = as.character(toKeep$Taxid))
allBac_clade	<- extract.clade(finalBacil_mp, node = bac_mrca)

# Then remove all the unnecessary internal tips
bacilPrune_data	<- drop.tip(allBac_clade, tip = as.character(toRemove$Taxid))
message(paste0("Group is monophyletic in Consensus: ", is.monophyletic(bacilPrune_data, tip = as.character(toKeep$Taxid))))
cons_mrca	<- getMRCA(bacilPrune_data, tip = as.character(toKeep$Taxid))
cons_clade	<- extract.clade(bacilPrune_data, node = cons_mrca)

# Check whether the ASTRAL phylogeny is correct
message(paste0("The astral phylogeny and consensus phylogeny for this group are identical: ",all.equal.phylo(astr_clade, cons_clade, use.edge.length = FALSE)))

# Finally write out the ASTRAL consensus tree (without the removed tips)
write.tree(astralPrune_data, file = file.path(newAstra_dir, "Astral_speciesTree_ultraChecked.txt"))