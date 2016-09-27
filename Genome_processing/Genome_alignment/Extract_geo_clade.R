### Extract the geobacillus clade from the Astral (final) species tree and write it out as seperate for potential seeding of the multi-genome alignment of Geobacillus ###
setwd("/users/aesin/Desktop/Mowgli/Species_tree")
spec_tree <- "Ultrametric_species_tree.txt"
spec_tree_phylo <- read.tree(spec_tree)

geobacs <- grep("Geobacillus",spec_tree_phylo$tip.label,value=T)
getMRCA(spec_tree_phylo, geobacs)
geo_extracted <- extract.clade(spec_tree_phylo, 5789)

write.tree(geo_extracted, file="/users/aesin/desktop/Geo_analysis/Geo_omes/Genome_alignment/Geobac_guide_tree/Species_tree_extracted_19_geo.txt")
## The re-read is used to convert spaces " " in tip labels into underscores "_" ##
re_read_geo_tree <- read.tree("/users/aesin/desktop/Geo_analysis/Geo_omes/Genome_alignment/Geobac_guide_tree/Species_tree_extracted_19_geo.txt")
geo_tips <- re_read_geo_tree$tip.label

write.table(geo_tips, file="Geo_tips.txt", sep="\n", quote=FALSE, col.names=FALSE, row.names=FALSE)

write.tree(geo_ext_no_edges, file="/users/aesin/desktop/Geo_analysis/Geo_omes/Genome_alignment/Geobac_guide_tree/Species_tree_extracted_19_geo_topo.txt")