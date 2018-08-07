#!/usr/bin/env Rscript
invisible(sapply(HGTPos.all, source, .GlobalEnv))

# Any libraries explicitly used in the script
require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load("tidyverse")

# General position data
perTypeData		<- readRDS(file.path(positionData_path, "AG_perTypeData.rds"))

# Open All_prot database
dbConn			<- dbConnect(RSQLite::SQLite(), allProtDB_path)

# -------------------------------------- #
# Get all AG prots from the database
allAG_prots		<- dbSendQuery(dbConn, 'SELECT * FROM t1 WHERE is_ag = :agFlag')
dbBind(allAG_prots, param = list(agFlag = 1))
allAG_prot_raw	<- dbFetch(allAG_prots)
dbClearResult(allAG_prots)

allAG_prot_df <- allAG_prot_raw
allAG_prot_df %<>% dplyr::rename(orthGroup = OrthGroup)
allAG_products <- allAG_prot_df %>% select(protID, product) %>% dplyr::rename(productFunction = product)
# %>% replace_na(list(orthGroup = "NA"))

coreColNames <- c("type", "protID", "orthGroup", "locusTag", "taxid", "binomial", "COGcat", "geneStart", "geneEnd", "strand", "productFunction", "relGeneStart", "relGeneEnd", "relStrand", "distToOri")

# -------------------------------------- #
# Define HGT genes
HGT_posEntries <- perTypeData$lHGT$`4`$allPosData
HGT_protIDs <- HGT_posEntries %>% pull(protID)

# Add "HGT Age", remove subgroup column, remove
HGT_final <- HGT_posEntries %>% 
	mutate(type = "HGT") %>%
	mutate("hgtAge" = case_when(
		Subgroup == TRUE ~ "Recent",
		TRUE ~ "Old")) %>%
	mutate(Subgroup = NULL, CircStart = NULL, CircEnd = NULL) %>%
	mutate(COGcat = unlist(lapply(COGcat, paste, collapse = "|"))) %>%
	dplyr::rename(receiverEdge = recepEdge) %>%
	left_join(allAG_products, by = "protID")

notCoreColNames <- names(HGT_final)[!names(HGT_final) %in% coreColNames]

HGT_final %<>% select(coreColNames, notCoreColNames)

# Remove the donor/receptor branch data
HGT_final %<>% select(-c(donorEdge, receiverEdge))

# Remove all HGT genes from overall set
allAG_prot_df %<>% filter(!protID %in% HGT_protIDs)

# -------------------------------------- #
# Define vertical genes
Ver_posEntries <- perTypeData$Ver$`3`$allPosData
Ver_protIDs <- Ver_posEntries %>% pull(protID)

Ver_final <- Ver_posEntries %>% 
	mutate(type = "Vertical") %>%
	mutate(CircStart = NULL, CircEnd = NULL) %>%
	mutate(COGcat = unlist(lapply(COGcat, paste, collapse = "|"))) %>%
	left_join(allAG_products, by = "protID") %>%
	select(coreColNames)

# Remove all HGT genes from overall set
allAG_prot_df %<>% filter(!protID %in% Ver_protIDs)

# -------------------------------------- #
# Define Grey Zone genes
all_posEntries <- perTypeData$All$allPosData
grey_posEntries <- all_posEntries %>% filter(!protID %in% Ver_protIDs & !protID %in% HGT_protIDs)
grey_protIDs <- grey_posEntries %>% pull(protID)

Grey_final <- grey_posEntries %>% 
	mutate(type = "Grey Zone") %>%
	mutate(CircStart = NULL, CircEnd = NULL) %>%
	mutate(COGcat = unlist(lapply(COGcat, paste, collapse = "|"))) %>%
	left_join(allAG_products, by = "protID") %>%
	select(coreColNames)


# Remove all HGT genes from overall set
allAG_prot_df %<>% filter(!protID %in% grey_protIDs)

# -------------------------------------- #
# Bind HGT, Vertical and Grey Zone
allChrom_posEntries <- bind_rows(list(HGT_final, Ver_final, Grey_final))

# -------------------------------------- #

AG_plasmid <- allAG_prot_df %>% 
	mutate(type = "Plasmid") %>%
	dplyr::rename(locusTag = locus, geneStart = gene_start, geneEnd = gene_end, productFunction = product) %>%
	select(type, protID, orthGroup, locusTag, taxid, binomial, COGcat, geneStart, geneEnd, strand, productFunction)

# -------------------------------------- #
# Bind all positional data with plasmid data
allAG_tableOut <- bind_rows(allChrom_posEntries, AG_plasmid)

allAG_tableOut_COGfix <- allAG_tableOut %>% mutate(COGcat = str_replace(COGcat, "-", "ND")) 
allAG_tableOut_COGfix <- allAG_tableOut_COGfix %>% 
	mutate_at(vars(COGcat), funs(ifelse(COGcat %in% c("", " ", "NA"), NA, COGcat)))

write.table(allAG_tableOut_COGfix, file = file.path(master_path, "..", "AG_manuscript", "Tables", "Sup_tbl_2.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)










