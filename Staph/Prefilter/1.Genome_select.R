#!/usr/bin/env Rscript

# Any libraries explicitly used in the script
require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load(tidyverse, fs, ape)

## ----------------------------------------------- ##

master			<- "/Users/aesin/Desktop/Staph"
genomeList_dir	<- file.path(master, "Genomes", "Genome_lists")
prefilter_dir	<- file.path(master, "Prefilter", "Staphylococcus")


accAss_taxid_tbl	<- read_tsv(file.path(genomeList_dir, "Acc_ass_taxid_table.tsv"), col_names = c("Acc_ass", "Taxid"))
taxidBinom_tbl		<- read_tsv(file.path(genomeList_dir, "Taxid_refinedBinomial_table.tsv"))

proteome_dir		<- file.path(master, "..", "Geo_again", "Proteomes")

fullAstral_tree		<- file.path(master, "..", "Geo_again", "Astral_speciesTree", "Correct_anogeo", "Astral_speciesTree_correct_ultra.txt")

## ----------------------------------------------- ##
# Write a binomial-tipped astral tree for visual analysis
fullAstral_data		<- read.tree(fullAstral_tree)

fullAstral_binom	<- fullAstral_data
fullAstral_binom$tip.label	<- left_join(tibble(Taxid = as.integer(fullAstral_binom$tip.label)), taxidBinom_tbl, by = "Taxid") %>% pull(Binomial)
write.tree(fullAstral_binom, file.path(master, "fullAstral_binomial.txt"))

# Extract the staphylococcaceae clade: Gemella and Staph aureus
edge_tips	<- c("1785995", "1280")
mrca_node	<- getMRCA(fullAstral_data, edge_tips)

get_clade	<- extract.clade(fullAstral_data, node = mrca_node)
staph_taxidBinom	<- tibble(Taxid = as.integer(get_clade$tip.label)) %>% left_join(taxidBinom_tbl, by = "Taxid")
staph_AccAss		<- staph_taxidBinom %>% left_join(accAss_taxid_tbl, by = "Taxid")


## ----------------------------------------------- ##
# Create output master dir
output_dir	<- file.path(prefilter_dir, "Proteomes")
dir_create(output_dir)

write_delim(staph_AccAss, path = file.path(prefilter_dir, "Acc_ass_list.txt"))

# Transfer the relevant proteomes (clean and dup)
clean_dir	<- file.path(output_dir, "Clean")
dup_dir		<- file.path(output_dir, "Dup")
dir_create(c(clean_dir, dup_dir))

invisible(lapply(staph_AccAss$Acc_ass, function(acc_ass) {

	prot_name	<- paste0(acc_ass, "_proteome.fasta.gz")
	clean_prot	<- file.path(proteome_dir, "Proteome_clean_fastas", prot_name)
	dup_prot	<- file.path(proteome_dir, "Proteome_dup_fastas", prot_name)

	file.copy(clean_prot, clean_dir)
	file.copy(dup_prot, dup_dir)
}))
