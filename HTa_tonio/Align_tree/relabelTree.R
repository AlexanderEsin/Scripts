#!/usr/bin/env Rscript

require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load("RSQLite", "tidyverse", "parallel", "seqinr", "ape", "magrittr")

master_dir		<- "/Users/aesin/Desktop/HTa_tonio/OLD/"
tree_dir		<- file.path(master_dir, "Alignment", "Raxml")

# Open All_prot database
HTa_db_path		<- file.path(master_dir, "HTa_prot_db")
dbConn			<- dbConnect(RSQLite::SQLite(), HTa_db_path)

## Read in nodes/names
taxdmp_path	<- file.path("/Users/aesin/Desktop/Geo_again", "Genomes", "taxdmp")
nodes_data	<- as_tibble(read_rds(file.path(taxdmp_path, "nodes.rds")))
names_data	<- as_tibble(read_rds(file.path(taxdmp_path, "names.rds")))

# Read in tree
tree_data		<- read.tree(file.path(tree_dir, "RAxML_bipartitions.raxml_tree.txt"))

# Get binomials for tree names
binomial_tbl	<- dbSendQuery(dbConn, 'SELECT protID, binomial, taxid FROM t1 WHERE protID = :protIDs')
dbBind(binomial_tbl, param = list(protIDs = tree_data$tip.label))
binomial_raw_df	<- dbFetch(binomial_tbl)
dbClearResult(binomial_tbl)

# New column with just first two words of the binomial name
binomial_df <- binomial_raw_df %>% select(-taxid) %>% mutate(binomial_trunc = word(binomial, 1, 2, sep = " "))
# Count the number of duplicate organism names, assign them a counter
binomial_df	%<>% group_by(binomial_trunc) %>% mutate(binom_count = seq(1:n()))
# Combine the binomial and counter so each label can be traced back to protID
binomial_df %<>% mutate(new_tip_label = paste(unlist(lapply(str_split(binomial_trunc, " "), paste, collapse = "_")), binom_count, sep = "_"))

write.table(binomial_df, file = file.path(tree_dir, "..", "protID_label_mapTbl.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

reLabelled_tree	<- tree_data
reLabelled_tree$tip.label	<- binomial_df$new_tip_label

write.tree(reLabelled_tree, file = file.path(tree_dir, "reLabelled_RAxML_HTa_tree.txt"))

dbDisconnect(dbConn)





unique_taxids	<- unique(binomial_raw_df$taxid)
domNames_l	<- lapply(1:length(unique_taxids), function(index) {
	
	
	taxid <- unique_taxids[index]
	message(paste0("Index: ", index, " Taxid: ", taxid))
	taxidNum <- as.numeric(taxid)

	while (TRUE) {

		if (identical(taxidNum, 1)) break

		entry	<- nodes_data %>% filter(taxid == taxidNum)
		rank	<- entry$nodeRank
		if (identical(rank, "superkingdom") || length(rank) == 0) {break}
		taxidNum	<- entry$parentTaxid
	}

	if (identical(taxidNum, 1) | length(rank) == 0) return(tibble(taxid = taxid, nodeRank = NA))

	return(tibble(taxid = taxid, nodeRank = taxidNum))
})

domNames_df <- bind_rows(domNames_l)

binomial_wRank		<- left_join(binomial_raw_df, domNames_df, by = "taxid")
# Manually change the 1 missing nodeRank - it's archaea
binomial_wRank[which(is.na(binomial_wRank$nodeRank)), "nodeRank"] <- 2157

binomial_wRankName	<- left_join(
	binomial_wRank,
	names_data %>%
		filter(nameClass == "scientific name")%>%
		select(taxid, name),
	by = c("nodeRank" = "taxid"))

binomial_wRankName <- as_tibble(binomial_wRankName)
# Check what we have
# binomial_wRankName %>% group_by(name) %>% summarise(n = n())
archIDs	<- binomial_wRankName %>% filter(name == "Archaea") %>% pull(protID)

eukIDs	<- binomial_wRankName %>% filter(name == "Eukaryota")



domNames_l	<- domNames_l[!is.na(domNames_l)]
if (length(domNames_l) == 0) {
	topFam <- NA
} else {
	famName_tbl	<- sort(table(unlist(famNames_l)), decreasing = T)
	topFam		<- names(famName_tbl)[1]
}






























