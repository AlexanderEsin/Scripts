#!/usr/local/bin/Rscript

library(ape)
library(RSQLite)
library(stringr)

## Command line arguments
args <- commandArgs(TRUE)
evalue			<- "1e-50"
trunc_eval		<- str_sub(evalue, 3, -1)

## Directory structures
direct			<- file.path("/Users/aesin/Desktop/Bacillus/Prefilter/Bacillus_withOut")
# raxml_tree_dir	<- file.path(direct, "Family_groups", trunc_eval, "RAxML_tree")
ft_tree_dir		<- file.path(direct, "Family_groups", trunc_eval, "FT_tree")

out_tax_dir		<- file.path(direct, "Final_trees/Taxid_labelled")
out_bin_dir		<- file.path(direct, "Final_trees/Binomial_labelled")

dir.create(out_tax_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_bin_dir, recursive = TRUE, showWarnings = FALSE)


## Open All_prot_db database
database_path	<- "/Users/aesin/Desktop/Bacillus/All_prot_db_new"
conn			<- dbConnect(RSQLite::SQLite(), database_path)

## Read in the reconstructed group tree which has taxid ID as labels
# rax_tree_file	<- file.path(raxml_tree_dir, paste0("RAxML_bipartitions.super_tree", trunc_eval, ".txt"))
ft_tree_file	<- file.path(ft_tree_dir, paste0("FT_super_tree", trunc_eval, ".txt"))
tree_list		<- list(FT = ft_tree_file)

tree_adj_tips_list	<- lapply(tree_list, function(group_tree_file) {
	group_tree_data	<- read.tree(group_tree_file)
	taxid_labels	<- group_tree_data$tip.label

	## For each taxid table, get the binomial and strain name from the database
	binomial_tbl	<- dbSendQuery(conn, 'SELECT binomial, strain FROM t1 WHERE taxid = :taxids LIMIT 1')
	dbBind(binomial_tbl, param = list(taxids = taxid_labels))
	name_df			<- dbFetch(binomial_tbl)
	dbClearResult(binomial_tbl)

	# All binomials are unique, no need to use strain names
	binomial_list	<- name_df$binomial
	if (length(unique(binomial_list)) != length(binomial_list)) {
		stop("Binomial names are not unique: error A1")
	}

	refined_binoms	<- lapply(binomial_list, function(name) {
		no_weird	<- gsub("[^A-Za-z0-9.[:space:]]", "", name)
		no_dbl_spc	<- gsub("[[:space:]]{2}", " ", no_weird)
		underscore	<- gsub("[[:space:]]", "_", no_dbl_spc)
		return(underscore)
	})

	# Still 5070 unique names
	if (length(unique(refined_binoms)) != length(refined_binoms)) {
		stop("Binomial names are not unique after name refinement: error A2")
	}

	file.copy(group_tree_file, out_tax_dir, overwrite = TRUE)

	## Write out the new binomial-tipped tree
	tip_change_tree				<- group_tree_data
	tip_change_tree$tip.label	<- refined_binoms
	return(tip_change_tree)
})
names(tree_adj_tips_list) <- names(tree_list)

for (i in 1:length(tree_adj_tips_list)) {
	tree_type	<- names(tree_adj_tips_list)[i]
	tree_data	<- tree_adj_tips_list[[i]]
	write.tree(tree_data, file = file.path(out_bin_dir, paste0(tree_type, "_binomial_SupTree_", trunc_eval, ".txt")))
}

dbDisconnect(conn)

