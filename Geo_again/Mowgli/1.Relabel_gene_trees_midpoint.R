#!/usr/local/bin/Rscript

## Relabel the FastTree gene trees for input to Mowgli:

# 1.	The species tree has 5070/5073 taxa in dataset. Need
#		to remove any occurence of the three taxa from all
#		gene

# 2.	We need to relabel all tips in the gene trees to just
#		the taxid - to match the species tree

# 3.	Any duplicate taxids in the gene tree tips (two+ proteins
#		from the same genome) must be labelled accordingly with the
#		"taxid_[num]" style for Mowgli and a key must be made between
#		full protID labels and the taxid style

## Packages
library(ape)
library(dplyr)
library(geiger)
library(gtools)
library(phangorn)
library(phytools)
library(RSQLite)
library(stringr)

## Paths and key files
direct			<- "/Users/aesin/Desktop/Geo_again"
speciesT_file	<- file.path(direct, "Astral_speciesTree", "Astral_noDup_speciesTree_ultra.txt")
database_path	<- file.path(direct, "All_prot_db")
geneT_input_dir	<- file.path(direct, "Group_fastTree", "Final_trees")
mowgli_out_dir	<- file.path(direct, "Mowgli", "GeneTree_input")

## Make output dir
dir.create(mowgli_out_dir, showWarnings = FALSE, recursive = TRUE)

## Open database connection
dbConn	<- dbConnect(RSQLite::SQLite(), database_path)

## Open species tree and get a list of unique tips
species_data	<- read.tree(speciesT_file)
species_taxids	<- species_data$tip.label
unique_taxids	<- length(unique(species_taxids))

## Get a list of the gene tree files
gene_tree_list	<- dir(geneT_input_dir, pattern = "*.txt", full.names = TRUE)
gene_tree_list	<- mixedsort(gene_tree_list)
# [6000:length(gene_tree_list)]
gene_tree_path	<- gene_tree_list[6446]

## Set up the counter
done_counter	<- 1

process_table_list	<- lapply(gene_tree_list, function(gene_tree_path) {
	## Process the file name to get the group number
	gene_tree_file	<- basename(gene_tree_path)
	group_number	<- str_sub(gene_tree_file, start = 1, end = -13)

	## Tracker
	message(paste0("Gene tree processed: ", group_number))

	## Read in gene tree and get list of tips
	gene_tree_data	<- read.tree(gene_tree_path)
	gene_tree_tips	<- gene_tree_data$tip.label
	original_tips	<- length(gene_tree_tips)

	## Translate the protID tips into corresponding taxids
	gene_taxid_tbl	<- dbSendQuery(dbConn, 'SELECT DISTINCT protID, taxid FROM t1 WHERE protID = :tips')
	dbBind(gene_taxid_tbl, param = list(tips = gene_tree_tips))
	gene_taxid_df	<- dbFetch(gene_taxid_tbl)
	dbClearResult(gene_taxid_tbl)

	## 1.	Make sure that the tree does not contain any of the
	##		3 taxa not represented in the species tree
	gene_tree_taxid	<- gene_taxid_df$taxid
	extra_taxids	<- gene_tree_taxid[which(!gene_tree_taxid %in% species_taxids)]
	if (length(extra_taxids) != 0) {
		# unique_extra_taxids	<- length(unique(extra_taxids))
		message(paste0("\t--> Unrepresented tips in gene tree: ", length(extra_taxids)))
		extra_protIDs	<- gene_taxid_df$protID[which(gene_taxid_df$taxid %in% extra_taxids)]
		
		# New tree data - with extra tips pruned
		gene_tree_data	<- drop.tip(gene_tree_data, tip = extra_protIDs)

		# Translation dataframe with extra tips removed
		gene_taxid_df	<- gene_taxid_df[which(!gene_taxid_df$taxid %in% extra_taxids),]
	}

	## Check the number tips on the tree equals number of rows in the translation table
	if (length(gene_tree_data$tip.label) != nrow(gene_taxid_df)) {
		stop(paste0("Number of tips doesn't equal the number of rows in translation table. Group: ", group_number))
	}
	## Check all protIDs are unique
	if (length(gene_taxid_df$protID) != length(unique(gene_taxid_df$protID))) {
		stop(paste0("Non unique protIDs present in group: ", group_number))
	}

	extra_pruned_tips	<- length(gene_tree_data$tip.label)

	## 2. Relabel the taxid labels so that multiple taxids belonging to the same species
	## are labelled accordingly

	## Each row correspond to one tip - relabel tip by tip. Check for multiples
	## at each row, and relabel them if necessary (*_1, *_2, *_3 ...)
	used_taxids	<- list()
	for (index in 1:nrow(gene_taxid_df)) {
		# Get row and protID/taxid
		row		<- gene_taxid_df[index,]
		protID	<- row$protID
		taxid	<- row$taxid

		# Check for any identical taxids already relabelled
		# Label next taxid accordingly
		test_number	<- 1
		test_taxid	<- paste0(taxid, "_", test_number)
		while(length(used_taxids[which(used_taxids == test_taxid)]) != 0) {
			test_number	<- test_number + 1
			test_taxid	<- paste0(taxid, "_", test_number)
		}

		new_taxid	<- test_taxid
		used_taxids	<- c(used_taxids, new_taxid)

		# Change the taxid in the table
		gene_taxid_df$taxid[which(gene_taxid_df$protID == protID)]	<- new_taxid
	}

	## Check there are no taxid duplicates (would lead to duplicate tips)
	if (nrow(gene_taxid_df) != length(unique(gene_taxid_df$taxid))) {
		stop("Non-unique taxids are present after exchange")
	}

	## Relabel tree
	gene_tree_data$tip.label	<- gene_taxid_df$taxid

	## Force binary, and midpoint
	## In some cases - e.g. 6148_FT_tree.txt - all the edges are 0
	## presumably alignment is identical. So, we can't midpoint and this tree cannot be used.
	gene_tree_di	<- multi2di(gene_tree_data)
	# If we had to resolve a polychotomy, then binary_resolved = TRUE
	binary_resolved	<- !all.equal.phylo(gene_tree_data, gene_tree_di)
	if (sum(node.depth.edgelength(gene_tree_di)) == 0) {
		message(paste0("\t--> All edges same length, skipping: ", group_number))
		return(data.frame(Group = group_number, Start.Tips = original_tips, Extra.Pruned = extra_pruned_tips, Resolved.Poly = binary_resolved, All.Equal.Branches = TRUE, Root.Algorithm = NA, stringsAsFactors = FALSE))
	}
	
	## Root the tree - choosing the correct algorithim (the fast phytools midpoint.root sometimes fails)
	gene_tree_mp	<- try(midpoint.root(gene_tree_di), silent = TRUE)
	if (class(gene_tree_mp) == "try-error") {
		gene_tree_mp	<- midpoint(gene_tree_di)
		root_alg		<- "phangorn" 
	} else {
		root_alg		<- "phytools"
	}

	## Convert the numeric (non-root) node labels (0-1 SH test) to 0-100 (BS-like)
	num_nodeLab_indx	<- grep("[^A-Za-z]",gene_tree_mp$node.label)
	gene_tree_mp$node.label[num_nodeLab_indx]	<- as.character(round(as.numeric(gene_tree_mp$node.label[num_nodeLab_indx]) * 100))

	## Write output files: the midpoint-rooted tree and key
	output_folder	<- file.path(mowgli_out_dir, group_number)
	dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

	write.tree(gene_tree_mp, file = file.path(output_folder, paste0(group_number, "_FT_relabelled.txt")))
	write.table(gene_taxid_df, file = file.path(output_folder, paste0(group_number, "_KEY_tips.txt")), row.names = FALSE, quote = FALSE, sep = "\t")

	return(data.frame(Group = group_number, Start.Tips = original_tips, Extra.Pruned = extra_pruned_tips, Resolved.Poly = binary_resolved, All.Equal.Branches = FALSE, Root.Algorithm = root_alg, stringsAsFactors = FALSE))
})

process_table_df	<- bind_rows(process_table_list)
write.table(process_table_df, file = file.path(direct, "Mowgli", "Preprocess_trees_stats.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)


dbDisconnect(dbConn)