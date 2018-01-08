#!/usr/local/bin/Rscript

library(ape)
library(RSQLite)

direct			<- "/Users/aesin/Desktop/Geo_again/Astral_speciesTree/"
species_tree	<- file.path(direct, "Astral_noDup_speciesTree_ultra.txt")
output_tree		<- file.path(direct, "Astral_noDup_refinedNames.txt")

## Open All_prot_db database
database_path	<- "/Users/aesin/Desktop/Geo_again/All_prot_db"
conn			<- dbConnect(RSQLite::SQLite(), database_path)

spec_tree_data	<- read.tree(species_tree)
taxid_labels	<- spec_tree_data$tip.label

## For each taxid table, get the binomial and strain name from the database
binomial_tbl	<- dbSendQuery(conn, 'SELECT binomial, strain FROM t1 WHERE taxid = :taxids LIMIT 1')
dbBind(binomial_tbl, param = list(taxids = taxid_labels))
name_df			<- dbFetch(binomial_tbl)
dbClearResult(binomial_tbl)

# All 5070 binomials are unique, no need to use strain names
binomial_list	<- name_df$binomial
length(unique(binomial_list))

refined_binoms	<- lapply(binomial_list, function(name) {
	no_weird	<- gsub("[^A-Za-z0-9.[:space:]]", "", name)
	no_dbl_spc	<- gsub("[[:space:]]{2}", " ", no_weird)
	underscore	<- gsub("[[:space:]]", "_", no_dbl_spc)
	return(underscore)
})

# Still 5070 unique names
length(unique(refined_binoms))

tip_change_tree	<- spec_tree_data
tip_change_tree$tip.label	<- refined_binoms

write.tree(tip_change_tree, file = output_tree)

## Write a table linking the refined binomial to the taxid
taxid_binom_df <- cbind(Taxid = taxid_labels, Binomial = refined_binoms)
write.table(taxid_binom_df, file = file.path("/Users/aesin/Desktop/Geo_again/Genomes/Genome_lists", "Taxid_refinedBinomial_table.tsv"), sep = "\t", quote = FALSE, row.name = FALSE)


dbDisconnect(conn)

