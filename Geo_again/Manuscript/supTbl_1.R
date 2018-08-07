#!/usr/bin/env Rscript

# Load master variables and HGT position functions
invisible(sapply(HGTPos.all, source, .GlobalEnv))

# Any libraries explicitly used in the script
require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load("RSQLite", "dplyr","ggplot2", "ggpubr", "wesanderson", "gridExtra", "parallel", "data.table")


# Open All_prot database
dbConn			<- dbConnect(RSQLite::SQLite(), allProtDB_path)

# ------------------------------------------------ #

# Open subgroup table
speciesSubgroups_file	<- file.path(genome_path, "Genome_lists", "Species_groupings.tsv")
speciesSubgroups_df		<- read.table(speciesSubgroups_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "")

# Get the Sqlite database data per taxon (1 return per species)
taxon_tbl	<- dbSendQuery(dbConn, 'SELECT * FROM t1 WHERE taxid = :taxids LIMIT 1')
dbBind(taxon_tbl, param = list(taxids = speciesSubgroups_df$taxid))
taxon_df	<- dbFetch(taxon_tbl)
dbClearResult(taxon_tbl)


# ------------------------------------------------ #

# Supergroups to label table
sporeGroups		<- c("Bacillaceae", "Paenibacillaceae", "Alicyclobacillaceae", "Pasteuriaceae", "Sporolactobacillaceae", "Thermoactinomycetaceae")
nonSporeGroups	<- c("Listeriaceae", "Staphylococcaceae", "Lactobacillales")

# Check which groups have not been done
allGroups		<- c(sporeGroups, nonSporeGroups)
groupsToPrepare	<- lapply(allGroups, function(taxonGroup) {
	file_name	<- paste0(taxonGroup, ".rds")
	file_path	<- file.path(supergroupData_path, file_name)
	if (!file.exists(file_path)) {
		return(taxonGroup)
	} else {
		return(NA)
	}
})
groupsToPrepare	<- groupsToPrepare[!is.na(groupsToPrepare)]

# If we need to prepare groups, read in the Nodes and Names data tables
if (length(groupsToPrepare) > 0) {
	nodes_dt	<- readRDS(file.path(taxdmp_path, "nodes.rds"))
	names_dt	<- readRDS(file.path(taxdmp_path, "names.rds"))

	# Prepare and save the groups that haven't been done before
	invisible(lapply(groupsToPrepare, function(taxonGroup) {

		# Identify the parentSpecies taxids that belong to desired superGroups
		superGroup_taxids	<- taxidsBySupergroup(taxonGroup, speciesGroupings = speciesSubgroups_df, nodes_dt = nodes_dt, names_dt = names_dt)

		# Save the result
		outputName	<- file_name	<- paste0(taxonGroup, ".rds")
		saveRDS(object = superGroup_taxids, file = file.path(supergroupData_path, outputName))

	}))

	# Clean away the names and nodes data tables
	rm(list = c("names_dt", "nodes_dt")); gc()

}


superGroup_readList	<- lapply(allGroups, function(taxonGroup) {
	file_name	<- paste0(taxonGroup, ".rds")
	file_path	<- file.path(supergroupData_path, file_name)

	taxonData	<- readRDS(file_path)
	return(taxonData)
})
names(superGroup_readList)	<- allGroups

# ------------------------------------------------ #

superGroup_df	<- bind_rows(lapply(1:length(superGroup_readList), function(listIndex) {
	famName		<- names(superGroup_readList)[listIndex]
	famTaxids	<- superGroup_readList[[listIndex]]

	out_df		<- data.frame(taxid = famTaxids, group = rep(famName, length(famTaxids)), stringsAsFactors = FALSE)
	return(out_df)
}))

# ------------------------------------------------ #


GPA_taxa	<- taxon_df %>% filter(is_ag == 1) %>%
	select(c(taxid, binomial, strain, acc_ass)) %>%
	left_join(speciesSubgroups_df, by = "taxid") %>%
	mutate(group = "GPA")

allTaxa		<- taxon_df %>% select(c(taxid, binomial, strain, acc_ass)) %>%
	filter(!taxid %in% GPA_taxa$taxid) %>%
	left_join(speciesSubgroups_df, by = "taxid") %>% 
	left_join(superGroup_df, by = c("parentSpecies" = "taxid")) %>%
	arrange(group) %>%
	mutate_at(.vars = vars(group), .funs = funs(ifelse(is.na(.), "ND", .)))

# Genome length extract
genome_l_df	<- taxon_df %>% select(c(taxid, genome_l))


# Bind the the Acc_ass to the speciesSubgroups table
supTbl_1	<- bind_rows(list(GPA_taxa, allTaxa)) %>%
	left_join(genome_l_df, by = "taxid") %>%
	dplyr::rename(
		accessionAssembly = acc_ass,
		parentTaxid = parentSpecies,
		parentBinomial = parentName,
		genomeLength = genome_l
	) %>%
	select(taxid, binomial, strain, accessionAssembly, parentTaxid, parentBinomial, genomeLength, group)

write.table(supTbl_1, file = file.path(master_path, "..", "AG_manuscript", "Tables", "Sup_tbl_1.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)