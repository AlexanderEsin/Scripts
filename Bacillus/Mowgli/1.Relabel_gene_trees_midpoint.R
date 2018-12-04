#!/usr/local/bin/Rscript

#	Any duplicate taxids in the gene tree tips (two+ proteins
#	from the same genome) must be labelled accordingly with the
#	"taxid_[num]" style for Mowgli and a key must be made between
#	full protID labels and the taxid style

# Any libraries explicitly used in the script
require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load(tidyverse, fs, magrittr, ape, phytools, phangorn, parallel)

# ----------------------------------------------------- #
# Paths and key files
direct			<- "/Users/aesin/Desktop/Bacillus"

speciesT_file	<- file.path(direct, "Astral_speciesTree", "Astral_speciesTree_ultraChecked.txt")

geneT_input_dir	<- file.path(direct, "Group_fastTree", "Final_trees")
mowgliOut_dir	<- file.path(direct, "Mowgli", "GeneTree_input")

# Taxid-binomial transl table
taxidBinom_tbl	<- read_tsv(file.path(direct, "Bac_genomes", "Genome_lists", "core_TaxidBinomial_table.tsv"))

# Make output dir
dir_create(mowgliOut_dir)

# Open species tree and get a list of unique tips
species_data	<- read.tree(speciesT_file)
species_taxids	<- species_data$tip.label
unique_taxids	<- length(unique(species_taxids))

# Get a list of the gene tree files
tree_tbl	<- tibble(filePath = dir(geneT_input_dir, pattern = "*.txt", full.names = TRUE)) %>%
	mutate(baseFileName = basename(filePath)) %>%
	mutate(groupNum = str_sub(baseFileName, start = 1, end = -13))

# ----------------------------------------------------- #

processed_list	<- mclapply(1:nrow(tree_tbl), function(rowIndex) {

	# Read in gene tree and get list of tips
	dataRow	<- tree_tbl[rowIndex,]

	geneTree_data	<- read.tree(dataRow$filePath)
	tip_tbl	<- tibble(rawGeneTreeTips = geneTree_data$tip.label)

	# Extract taxid from the protID
	tip_tbl	%<>% mutate(Taxid = str_extract(rawGeneTreeTips, "(?<=.)[0-9]+(?=_)"))

	# Check that the gene tree does not contain any tips not present in species tree
	extra_taxids	<- tip_tbl %>% filter(!Taxid %in% species_taxids)
	if (nrow(extra_taxids) > 0) {
		stop(paste0("Do not expect extra taxids. Group = ", dataRow$groupNum))
	}

	# Check all protIDs are unique
	if (nrow(tip_tbl) != nrow(distinct(tip_tbl, rawGeneTreeTips))) {
		stop(paste0("Non unique protIDs present. Group = ", dataRow$groupNum))
	}

	# Relabel tips so that each taxid has a unique number ID - for multiples, this number increases
	tip_tbl %<>% group_by(Taxid) %>%
		mutate(taxidRelab = paste0(Taxid, "_", seq(1:length(Taxid)))) %>%
		ungroup()

	# Relabel tree tips
	geneTree_data$tip.label	<- tip_tbl$taxidRelab

	# Resolve polychotomies if necessary
	geneTreeDi_data	<- multi2di(geneTree_data)

	# If we had to resolve a polychotomy, then binary_resolved = TRUE
	binary_resolved	<- !all.equal.phylo(geneTree_data, geneTreeDi_data)

	# If all the edge lengths are zero, tree is useless - return
	if (sum(node.depth.edgelength(geneTreeDi_data)) == 0) {
		out_tbl	<- tibble(Group = dataRow$groupNum, numTips = nrow(tip_tbl), resolvedPoly = binary_resolved, allEqualBranches = TRUE, rootAlgorithm = NA)
		return(out_tbl)
	}

	# Proceed to root tree - if the fast phytools midpoint.root fails, use phangorn
	geneTreeMP_data	<- try(midpoint.root(geneTreeDi_data), silent = TRUE)
	if (class(geneTreeMP_data) == "try-error") {
		geneTreeMP_data	<- midpoint(geneTreeDi_data)
		root_alg		<- "phangorn" 
	} else {
		root_alg		<- "phytools"
	}

	# Convert the numeric (non-root) node labels (0-1 SH test) to 0-100 (BS-like)
	numericNodeLab_index	<- which(str_detect(geneTreeMP_data$node.label, "[^A-Za-z]"))
	geneTreeMP_data$node.label[numericNodeLab_index]	<- as.character(round(as.numeric(geneTreeMP_data$node.label[numericNodeLab_index]) * 100))

	# Create output directory
	out_dir	<- file.path(mowgliOut_dir, dataRow$groupNum)
	dir_create(out_dir)

	# Write output files: the midpoint-rooted tree and key
	write.tree(geneTreeMP_data, file = file.path(out_dir, paste0(dataRow$groupNum, "_FT_relabelled.txt")))
	write_tsv(tip_tbl, path = file.path(out_dir, paste0(dataRow$groupNum, "_KEY_tips.txt")))

	# Return outout table
	out_tbl	<- tibble(Group = dataRow$groupNum, numTips = nrow(tip_tbl), resolvedPoly = binary_resolved, allEqualBranches = FALSE, rootAlgorithm = root_alg)
	return(out_tbl)

}, mc.cores = 15)

# Write out the summary table
processed_tbl	<- bind_rows(processed_list)
write_tsv(processed_tbl, path = file.path(direct, "Mowgli", "Preprocess_trees_stats.tsv"))
