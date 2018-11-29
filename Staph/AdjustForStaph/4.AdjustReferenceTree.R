#!/usr/bin/env Rscript

# Any libraries explicitly used in the script
require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load(tidyverse, fs, ape, phangorn)

master			<- file.path("/Users/aesin/Desktop/Staph")
core_list_dir	<- file.path(master, "Staph_genomes", "Genome_lists")

# Full (final) Astral species tree from the AG set
fullAstral_tree	<- file.path(master, "..", "Geo_again", "Astral_speciesTree", "Correct_anogeo", "Astral_speciesTree_correct_ultra.txt")
fullAstral_data	<- read.tree(fullAstral_tree)

# New Astral reference tree for this group reconciliation
newAstra_dir	<- file.path(master, "Astral_speciesTree")
dir_create(newAstra_dir)

# The consensus tree
finalCore_tree	<- file.path(master, "Prefilter", "Staphylococcus", "Final_trees", "Taxid_labelled", "FT_super_tree-50.txt")
finalCore_data	<- read.tree(finalCore_tree)

# Tables to read in
toRemove			<- read_tsv(file.path(core_list_dir, "coreToRemove.tsv"))
toKeep				<- read_tsv(file.path(core_list_dir, "coreToKeep.tsv"))

# ----------------------------------------------------- #
# First remove the excess tips that are no longer included in the dataset (extra core genomes)
astralPrune_data	<- drop.tip(fullAstral_data, tip = as.character(toRemove$Taxid))

# Check whether our core is monophyletic & extract clude from Astral tree
message(paste0("Group is monophyletic in Astral: ", is.monophyletic(astralPrune_data, tip = as.character(toKeep$Taxid))))
astr_mrca	<- getMRCA(astralPrune_data, tip = as.character(toKeep$Taxid))
astr_clade	<- extract.clade(astralPrune_data, node = astr_mrca)

# Now extract the relevant tips from the Core consensus tree
# First find the root node of our core group (the consensus tree has na outgroup)
finalCore_mp	<- midpoint(finalCore_data)
core_mrca		<- getMRCA(finalCore_mp, tip = as.character(toKeep$Taxid))
allCore_clade	<- extract.clade(finalCore_mp, node = core_mrca)

# Then remove all the unnecessary internal tips
corePrune_data	<- drop.tip(allCore_clade, tip = as.character(toRemove$Taxid))
message(paste0("Group is monophyletic in Consensus: ", is.monophyletic(corePrune_data, tip = as.character(toKeep$Taxid))))
cons_mrca	<- getMRCA(corePrune_data, tip = as.character(toKeep$Taxid))
cons_clade	<- extract.clade(corePrune_data, node = cons_mrca)

# Check whether the ASTRAL phylogeny is correct
message(paste0("The astral phylogeny and consensus phylogeny for this group are identical: ", all.equal.phylo(astr_clade, cons_clade, use.edge.length = FALSE)))

# There is a mismatch between the two phylogenies - in the consensus tree Taxids 61015 and 214473 (S. succinus and S. equorum) cluster together.
# In the astral phylogeny they are both independent groups.
# Here fix the Astral phylogeny.
consSubclade_fix	<- extract.clade(cons_clade, node = getMRCA(cons_clade, tip = as.character(c(61015, 1288))))
prune_targets		<- consSubclade_fix$tip.label
temp_names			<- paste0("BIND", seq(1:length(prune_targets)))

astralTest			<- astralPrune_data
# almost_pruned		<- drop.tip(astralPrune_data, anogeo_subset_phylo$tip.label[-1:-2])
astralTest$tip.label[which(astralTest$tip.label %in% prune_targets)] <- temp_names

rebind_where		<- getMRCA(astralTest, tip = temp_names)
rebind_subclade		<- bind.tree(astralTest, consSubclade_fix, where = rebind_where, position = 0.1)
rebind_subclade		<- multi2di(rebind_subclade, random = FALSE)
rebind_subclade$node.label	<- NULL

rebind_prune_keys	<- drop.tip(rebind_subclade, temp_names)

testFix_mrca	<- getMRCA(rebind_prune_keys, tip = as.character(toKeep$Taxid))
testFix_clade	<- extract.clade(rebind_prune_keys, node = testFix_mrca)

message(paste0("The fixed astral phylogeny and consensus phylogeny for this group are identical: ", all.equal.phylo(testFix_clade, cons_clade, use.edge.length = FALSE)))

# Finally write out the ASTRAL consensus tree (without the removed tips). This needs to be made into a cladogram with FigTree
write.tree(rebind_prune_keys, file = file.path(newAstra_dir, "Astral_speciesTree_Checked.txt"))

















