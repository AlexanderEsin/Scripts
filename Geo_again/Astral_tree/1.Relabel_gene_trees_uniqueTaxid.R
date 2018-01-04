#!/usr/local/bin/Rscript

library(ape)
library(RSQLite)
library(stringr)

direct			<- "/Users/aesin/Desktop/Geo_again/Group_fastTree/"
in_tree_dir		<- file.path(direct, "Final_trees")
out_tree_dir	<- file.path(direct, "Final_trees_relab")
no_dup_tree_dir	<- file.path(direct, "Final_trees_noDup")

database_path	<- "/Users/aesin/Desktop/Geo_again/All_prot_db"
dir.create(out_tree_dir, showWarnings = FALSE)
dir.create(no_dup_tree_dir, showWarnings = FALSE)

input_tree_list	<- dir(in_tree_dir, pattern = "*.txt", full.names = TRUE)

## Connect to database
conn	<- dbConnect(RSQLite::SQLite(), database_path)
message("Calculating number of unique taxids in the dataset...")
total_taxids		<- dbGetQuery(conn,'SELECT DISTINCT taxid from t1')
total_taxids		<- total_taxids$taxid
num_total_taxids	<- length(total_taxids)

## Number of taxa represented
all_taxids_in_trees		<- list()
all_taxids_represented	<- list()

done_counter	<- 1
for (input_tree in input_tree_list) {

	## Read in tree and get tips
	tree_data	<- read.tree(input_tree)
	tree_tips	<- tree_data$tip.label

	## In some circumastances, there are multiple entries for a particular protID in the database
	## This is because some genomes contain multiple plasmids each with an identical gene / protein at the same site
	## We will filter those out downstream, for now - limit to just 1

	## Query the Sqlite database using the tips as protID query
	taxid_tbl	<- dbSendQuery(conn, 'SELECT DISTINCT protID, taxid FROM t1 WHERE protID = :tips')
	dbBind(taxid_tbl, param = list(tips = tree_tips))
	taxid_df	<- dbFetch(taxid_tbl)
	dbClearResult(taxid_tbl)

	## The db results is a DF of a single column
	taxid_list	<- taxid_df$taxid
	if (length(taxid_list) != length(tree_tips)) {
		stop(paste0("The number of taxids retrieved does not equal the number of tips input for tree: ", input_tree))
	}
	
	## Produce the new tree object with taxids as tip labels
	tree_relab	<- tree_data
	tree_relab$tip.label	<- taxid_list

	## Define a new name for the new tree and write it out
	tree_basename	<- basename(input_tree)
	group_number	<- str_sub(tree_basename, start = 1, end = -13)
	new_file_name	<- paste0(group_number, "_FT_relab_tree.txt")
	write.tree(tree_relab, file = file.path(out_tree_dir, new_file_name))

	## Add a list of taxids represented in trees
	unique_taxids	<- unique(taxid_list)
	if (length(all_taxids_in_trees) != num_total_taxids) {
		all_taxids_in_trees	<- c(all_taxids_in_trees, unique_taxids)
		all_taxids_in_trees	<- unique(all_taxids_in_trees)
	}
	

	## Check whether all the taxids are unique in tree (i.e. no paralogs at all)
	num_unique_taxids	<- length(unique(taxid_list))
	if (length(taxid_list) == num_unique_taxids && length(tree_tips) >= 50) {
		write.tree(tree_relab, file = file.path(no_dup_tree_dir, new_file_name))

		all_taxids_represented	<- c(all_taxids_represented, unique_taxids)
		all_taxids_represented	<- unique(all_taxids_represented)
	}

	message(paste0("Relabelled: ", done_counter, " // ", length(input_tree_list), "... Taxids in all trees: ", length(all_taxids_in_trees), " // ", num_total_taxids, "... Taxids in noDup trees: ", length(all_taxids_represented), " // ", num_total_taxids))
	done_counter	= done_counter + 1
}

missing_taxids	<- setdiff(total_taxids, all_taxids_represented)
message(paste0("Taxids missing from the noDup tree set: ", paste(missing_taxids, collapse = " | ")))
missing_taxid_data	<- dbGetQuery(conn, 'SELECT acc_ass, binomial FROM t1 WHERE taxid = :taxids LIMIT 1', params = list(taxids = missing_taxids))
message(paste0("These correspond to: ", paste(missing_taxid_data$acc_ass, collapse = " | ")))
message(paste0("These correspond to: ", paste(missing_taxid_data$binomial, collapse = " | ")))

dbDisconnect(conn)
