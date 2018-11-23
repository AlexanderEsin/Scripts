#!/usr/bin/env Rscript

library(pacman)
p_load(tidyverse, seqinr, magrittr, parallel, dendextend)

# Process fasta header
processFastaHeader <- function(fastaVector) {

	byID_split		<- as_tibble(str_split_fixed(fastaVector, " ", 2)) %>%
		dplyr::rename(ID = 1) %>%
		mutate(ID = str_replace(ID, ">", ""))

	byName_split	<- byID_split %>%
		separate(2, into = c("Function", "Species"), sep = "\\[") %>%
		mutate(Species = str_replace(Species, "\\]", "")) %>%
		mutate(Function = str_trim(Function))

	processMulti	<- byName_split %>%
		mutate(Multispecies = case_when(
			str_detect(Function, "MULTISPECIES") ~ TRUE,
			TRUE ~ FALSE))

	processMulti	<- processMulti %>%
		mutate(Function = str_trim(str_replace(Function, "MULTISPECIES:", "")))

	return(processMulti)
}

# Assign taxid
assignTaxid	<- function(headerDF, names_data) {

	withTaxid	<- left_join(headerDF, names_data, by = c("Species" = "name")) %>%
		select(-uniqueName, -nameClass)

	# For unidentifed taxids, take first name of the species label (genus or higher classification) and find taxid
	notFound	<- withTaxid %>%
		filter(is.na(taxid)) %>%
		select(-taxid) %>%
		mutate(highClass = str_split_fixed(Species, " ", 2)[,1]) %>%
		left_join(., names_data, by = c("highClass" = "name")) %>%
		select(ID, taxid)

	# Where there are taxid NAs in the original table replace with taxid from the "notFound" table
	withTaxid %<>%
		left_join(notFound, by = "ID") %>%
		mutate(taxid = case_when(
			is.na(taxid.x) ~ taxid.y,
			TRUE ~ taxid.x)) %>%
		select(-taxid.x, -taxid.y)

	return(withTaxid)
}

# Process domains
processDomains	<- function(domains_tbl, headerDF) {

	domains_fix <- domains_tbl %>% 
		mutate(Query = str_trim(str_extract(Query, ">+.+?\\s+?"))) %>%
		mutate(Query = str_replace(Query, ">", "")) %>%
		dplyr::rename(ID = Query) %>%
		unite(combPos, From, To, Incomplete, sep = " ", remove = FALSE) %>%
		mutate(combPos = paste0("[", combPos, "]"))

	# For the combined dataframe
	by_nDomain	<- domains_fix %>%
		group_by(ID) %>%
		summarise(numDoms = n(), domPos = paste(combPos, collapse = ";"))
	combined_df	<- left_join(headerDF, by_nDomain, by = "ID")

	if (max(by_nDomain$numDoms) > 2) stop("There are more than 2 fold domains here .. code will not work")

	# For the split data frame
	domains_trim	<- domains_fix %>%
		select(ID, combPos, From, To, Incomplete) %>%
		left_join(headerDF, by = "ID") %>%
		group_by(ID) %>%
		arrange(From) %>%
		mutate(Section = case_when(
			n() == 1 ~ "Full",
			From == min(From) ~ "Nterm",
			TRUE ~ "Cterm")) %>%
		mutate(newID = paste(ID, Section, sep = "_")) %>%
		arrange(ID)

	return(list(combPos = combined_df, splitPos = domains_trim))
}

# Split fasta sequences by domain
splitFastaByDomain <- function(fastaSeqList, headerDF) {

	splitFasta_list	<- vector("list", nrow(headerDF))
	i	<- 1

	for (seqIndex in 1:length(fastaSeqList)) {

		seqName		<- names(fastaSeqList)[seqIndex]
		seq			<- fastaSeqList[[seqIndex]]
		attrName	<- attributes(seq)$name

		if (!identical(seqName, attrName)) stop(paste("ERROR: object and attribute names are not the same. Seq = ", seq))

		sections	<- headerDF %>% filter(ID == seqName)

		if (nrow(sections) == 2) {
			# Define the splitting position as the midpoint between the end of the first fold and the start of the second
			split_pos	<- round(min(sections$To) + ((max(sections$From) - min(sections$To)) / 2), 1)

			# Split the sequence by position
			nterm	<- seq[1:split_pos]
			cterm	<- seq[(split_pos + 1):length(seq)]

			# Cterm / Nterm names
			nterm_name	<- sections %>% filter(Section == "Nterm") %>% pull(newID)
			cterm_name	<- sections %>% filter(Section == "Cterm") %>% pull(newID)

			# Assign the nterm / cterm sequences to their own SeqFastaAA objects
			nterm_seq	<- as.SeqFastaAA(nterm, name = nterm_name)
			cterm_seq	<- as.SeqFastaAA(cterm, name = cterm_name)

			# Assign sequences to the list
			splitFasta_list[[i]]	<- nterm_seq
			splitFasta_list[[i+1]]	<- cterm_seq

			# Name list elements after sequences
			names(splitFasta_list)[i:(i+1)] <- c(nterm_name, cterm_name)

			i <- i + 2	
		} else if (nrow(sections) == 1) {
			# Get new ID from table
			newName	<- sections %>% pull(newID)

			# Set new ID as the name of the sequence (attr)
			attributes(seq)$name	<- newName

			# Add it to list
			splitFasta_list[[i]]	<- seq

			# Name the list element
			names(splitFasta_list)[i]	<- newName
			i	<- i + 1
		} else {
			stop("Some weird number of rows...")
		}
	}

	return(splitFasta_list)
}

# Reformat the header for the reformatted fasta
reformatFastaHeader <- function(headerDF) {
	# Reformat, and unite all with pipes
	heading <- headerDF %>% 
		mutate(Multispecies = case_when(
			Multispecies ~ "MULTISPECIES",
			TRUE ~ "SINGLE")) %>%
		unite(Heading, ID, Multispecies, Function, Species, taxid, numDoms, domPos, sep = " | ", remove = TRUE)

	withHeader <- bind_cols(list(headerDF, heading))

	return(withHeader)
}

# Reformat header for fasta split by domain
reformatSplitFastaHeader <- function(headerDF) {
	# Reformat, and unite all with pipes
	heading <- headerDF %>% 
		mutate(Multispecies = case_when(
			Multispecies ~ "MULTISPECIES",
			TRUE ~ "SINGLE")) %>%
		unite(Heading, newID, Multispecies, Function, Species, taxid, combPos, sep = " | ", remove = TRUE)

	withHeader <- bind_cols(list(headerDF, heading))

	return(withHeader)
}

# Set the new header on the fasta sequences
setNewHeader <- function(fastaSeqList, headerDF, split = FALSE) {

	newSeqAnnot_list <- mclapply(1:length(fastaSeqList), function(seqIndex) {

		seqName		<- names(fastaSeqList)[seqIndex]
		seq			<- fastaSeqList[[seqIndex]]
		attrName	<- attributes(seq)$name

		if (!identical(seqName, attrName)) stop("ERROR: object and attribute names are not the same")

		# Get the relevant (updated) header, and assign it as the "Annot" attribute to the sequence
		if (split) {
			matchHeader <- headerDF %>% filter(newID == seqName) %>% pull(Heading)
		} else {
			matchHeader <- headerDF %>% filter(ID == seqName) %>% pull(Heading)
		}
		
		attributes(seq)$Annot	<- matchHeader

		return(seq)
	}, mc.cores = 10)

	newNames	<- lapply(newSeqAnnot_list, function(seq) {
		name	<- attributes(seq)$Annot
		return(name)
	})

	names(newSeqAnnot_list)	<- newNames
	return(newSeqAnnot_list)
}

# For archaea, viruses and eukaryotes
processSequenceSet	<- function(taxGroup, out_dir, names_dt) {

	fasta	<- read.fasta(file = file.path(raw_dir, paste0("HU_IHF_", taxGroup, ".fasta")), seqtype = "AA")
	doms	<- read_tsv(file = file.path(dom_dir, paste0("domains_", taxGroup, ".txt")), comment = "")

	annot	<- unlist(getAnnot(fasta))

	# Process the header
	head_proc	<- processFastaHeader(annot)

	# Remove any sequences which are not formatted in the standard NCBI way
	toRemove	<- which(is.na(head_proc$Species))
	if (length(toRemove) != 0) {
		head_proc	<- head_proc[-toRemove,]
		fasta		<- fasta[-toRemove]
		doms		<- doms[-toRemove,]
	}
	

	# Get the taxids
	taxids	<- assignTaxid(head_proc, names_dt)
	# Process the domains data
	wDoms	<- processDomains(doms, taxids)

	# Write out the reformatted fasta - not split by domain (but containing domain info)
	header		<- reformatFastaHeader(wDoms$combPos)
	new_header	<- setNewHeader(fasta, header)
	write.fasta(new_header, names = names(new_header), file.out = file.path(out_dir, paste0("HU_IHF_", taxGroup, "_wTaxid.fasta")))

	# Split the sequences containing multiple domains
	splitByDom	<- splitFastaByDomain(fasta, wDoms$splitPos)
	splitHeader	<- reformatSplitFastaHeader(wDoms$splitPos)
	splitNewHd	<- setNewHeader(splitByDom, splitHeader, split = TRUE)

	# We will split into single domains (complete), single domains (incomplete), and double domains (all)
	singleDoms_IDs	<- wDoms$combPos %>% filter(numDoms == 1) %>% pull(ID)
	singleDom_cmpl	<- splitHeader %>% filter(ID %in% singleDoms_IDs & Incomplete == "-") %>% pull(Heading)
	singleDom_incm	<- splitHeader %>% filter(ID %in% singleDoms_IDs & Incomplete != "-") %>% pull(Heading)
	multiDom		<- splitHeader %>% filter(!ID %in% singleDoms_IDs) %>% pull(Heading)

	# Write the split fasta files
	singleDom_cmpl_fasta	<- splitNewHd[which(names(splitNewHd) %in% singleDom_cmpl)]
	singleDom_incm_fasta	<- splitNewHd[which(names(splitNewHd) %in% singleDom_incm)]
	multiDom_fasta			<- splitNewHd[which(names(splitNewHd) %in% multiDom)]

	write.fasta(splitNewHd, names = names(splitNewHd), file.out = file.path(out_dir, paste0("HU_IHF_", taxGroup, "_byDomSplit_all.fasta")))
	write.fasta(singleDom_cmpl_fasta, names = names(singleDom_cmpl_fasta), file.out = file.path(out_dir, paste0("HU_IHF_", taxGroup, "_oneDom_complete.fasta")))
	write.fasta(singleDom_incm_fasta, names = names(singleDom_incm_fasta), file.out = file.path(out_dir, paste0("HU_IHF_", taxGroup, "_oneDom_incomplete.fasta")))
	write.fasta(multiDom_fasta, names = names(multiDom_fasta), file.out = file.path(out_dir, paste0("HU_IHF_", taxGroup, "_twoDom_all.fasta")))

	out_tbl	<-  wDoms$splitPos %>% mutate(newTipName = paste(Species, newID, paste0("[", taxGroup, "]"), sep = " "))

	return(out_tbl)
}

# Once all bacteria are sorted by family, we can take up to 10 IDs from each one based on their pairwise alignment identity
# 1. Calculate pairwise identity matrix for all sequence in family, cluster them and select up to 10 clusters (based on separation)
# 2. Sample 1 sequence from each cluster
getFamSubset	<- function(famName) {
	bact_famSet	<- bact_allFamily %>% filter(familyName == famName)
	bact_famIDs	<- bact_famSet %>% pull(ID)
	bact_famSeq	<- bact_AAset[which(names(bact_AAset) %in% bact_famIDs)]

	numClusters	<- ifelse(nrow(bact_famSet) > 10, 10, nrow(bact_famSet))

	message(paste("\t-->", "Sequence number:", nrow(bact_famSet), "Clusters:", numClusters, sep = " "))

	if (nrow(bact_famSet) == 1) {
		fam_subset <- bact_famSet %>%
			mutate(Cluster = 1) %>%
			mutate(newName = paste(familyName, Cluster, sep = "_"))

		return(fam_subset)
	}

	dist		<- stringDist(bact_famSeq)
	clusters	<- dist %>%
		hclust(method = "complete") %>%
		dendextend::cutree(k = numClusters) %>%
		as.matrix(clusters) %>%
		as_tibble(rownames = "ID") %>%
		dplyr::rename(Cluster = 2)

	pickIDs		<- clusters %>% group_by(Cluster) %>% sample_n(1, replace = FALSE)

	fam_subset	<- bact_famSet %>%
		filter(ID %in% pickIDs$ID) %>%
		left_join(pickIDs, by = "ID") %>%
		mutate(newName = paste(familyName, Cluster, sep = "_"))

	return(fam_subset)
}



## ----------------------------------------------- ##

HTa_dir	<- "/Users/aesin/Desktop/HTa_tonio/Round2"
txd_dir	<- file.path(HTa_dir, "taxdmp")
seq_dir	<- file.path(HTa_dir, "Sequences")


raw_dir	<- file.path(seq_dir, "Raw")
dom_dir	<- file.path(seq_dir, "Domains")
out_dir	<- file.path(seq_dir, "Reformat")
tbl_dir	<- file.path(seq_dir, "Tables")

if (!dir.exists(out_dir)) dir.create(out_dir)
if (!dir.exists(tbl_dir)) dir.create(tbl_dir)

# Output fasta directories for each organism type
arch_dir	<- file.path(out_dir, "Archaea")
vir_dir		<- file.path(out_dir, "Viruses")
bact_dir	<- file.path(out_dir, "Bacteria")
euk_dir		<- file.path(out_dir, "Eukaryotes")

for (dir in c(arch_dir, vir_dir, bact_dir, euk_dir)) if (!dir.exists(dir)) dir.create(dir)



## ----------------------------------------------- ##
# Read taxdmp files

nodes_data	<- as_tibble(read_rds(file.path(txd_dir, "nodes.rds")))
names_data	<- as_tibble(read_rds(file.path(txd_dir, "names.rds")))
names_dataFilt <- names_data %>% group_by(name) %>% filter(row_number(taxid) <= 1) %>% ungroup()

## ----------------------------------------------- ##
# Archaea run
arch_tbl	<- processSequenceSet(taxGroup = "archaea", out_dir = arch_dir, names_dt = names_dataFilt)

# Virus run
vir_tbl	<- processSequenceSet(taxGroup = "viruses", out_dir = vir_dir, names_dt = names_dataFilt)

# Eukaryote run
euk_tbl	<- processSequenceSet(taxGroup = "euk", out_dir = euk_dir, names_dt = names_dataFilt)

# Write out all the non-Bacterial entries to a table to rename tree tips later
nonBact_tbl	<- bind_rows(list(arch_tbl, vir_tbl, euk_tbl))
write_tsv(nonBact_tbl, path = file.path(tbl_dir, "allNonBact_table.tsv"))

## ----------------------------------------------- ##


# Read bacteria fasta into AA string set
bact_AAset	<- readAAStringSet(file = file.path(raw_dir, "HU_IHF_bacteria.fasta"))
names(bact_AAset)	<- str_split_fixed(names(bact_AAset), " ", 2)[,1]

# Read bacteria fasta into seqinr list
bact_fasta	<- readAAStringSet(file = file.path(raw_dir, "HU_IHF_bacteria.fasta"))
bact_annot	<- unlist(getAnnot(bact_fasta))

# Process the headers
bact_proc	<- processFastaHeader(bact_annot)
# Remove any empty species column entries (~ 178 for bacteria)
bact_clean	<- bact_proc %>% filter(!str_length(Species) == 0 | is.na(Species))
bact_wTaxid	<- assignTaxid(bact_clean, names_dataFilt)
bact_wRank	<- bact_wTaxid %>% left_join(nodes_data %>% select(taxid, nodeRank), by = "taxid")


## Identify all bacteria to genus

# Extract just the sub-genus levels and identify the genus level.
bact_species	<- bact_wRank %>% filter(nodeRank %in% c("species", "no rank"))

genusData_df	<- lapply(1:nrow(bact_species), function(index) {

	print(index)

	row			<- bact_species[index,]
	taxidNum	<- row$taxid

	while (TRUE) {
		if (identical(taxidNum, 1)) break

		entry	<- nodes_data %>% filter(taxid == taxidNum)
		rank	<- entry$nodeRank
		if (identical(rank, "genus")) break
		taxidNum	<- entry$parentTaxid
	}

	if (identical(taxidNum, 1)) return(tibble(genusTaxid = NA, genusName = NA))

	genusTaxid	<- entry$taxid
	genusName	<- names_data %>% filter(taxid == genusTaxid & nameClass == "scientific name") %>% pull(name)

	return(tibble(genusTaxid = genusTaxid, genusName = genusName))
})

genusData_df	<- bind_rows(genusData_df)
bact_wGenus		<- bind_cols(bact_species, genusData_df)

# Add Genus column to those bacteria already at a genus-level
bact_genus	<- bact_wRank %>% filter(nodeRank == "genus") %>% mutate(genusTaxid = taxid, genusName = Species)
bact_allGenus	<- bind_rows(list(bact_genus, bact_wGenus))


familyData_df	<- lapply(1:nrow(bact_allGenus), function(index) {

	print(index)

	row			<- bact_allGenus[index,]
	taxidNum	<- row$genusTaxid
	if (is.na(taxidNum)) taxidNum <- row$taxid

	while (TRUE) {
		if (identical(taxidNum, 1)) break

		entry	<- nodes_data %>% filter(taxid == taxidNum)
		rank	<- entry$nodeRank
		if (identical(rank, "family")) break
		taxidNum	<- entry$parentTaxid
	}

	if (identical(taxidNum, 1)) return(tibble(familyTaxid = NA, familyName = NA))

	familyTaxid	<- entry$taxid
	familyName	<- names_data %>% filter(taxid == familyTaxid & nameClass == "scientific name") %>% pull(name)

	return(tibble(familyTaxid = familyTaxid, familyName = familyName))
})

# Add family column
familyData_df	<- bind_rows(familyData_df)
bact_wFamily	<- bind_cols(bact_allGenus, familyData_df)

# Add family column to those bacteria already at a family-level
bact_family		<- bact_wRank %>% filter(nodeRank == "family") %>% mutate(genusTaxid = NA, genusName = NA, familyTaxid = taxid)
bact_family		<- left_join(bact_family, names_data %>% filter(nameClass == "scientific name") %>% select(taxid, name), by = c("familyTaxid" = "taxid"))
bact_family		<- bact_family %>% dplyr::rename(familyName = name)

# All our bacterial proteins are classified to family (where possible)
bact_allFamily	<- bind_rows(list(bact_family, bact_wFamily))


# Remove an NA family names - that leaves 43,454 bact sequences
bact_allFamClean	<- bact_allFamily %>% filter(!is.na(familyName))
bact_famNames		<- bact_allFamClean %>% distinct(familyName) %>% pull(familyName)

write_tsv(bact_allFamClean, path = file.path(tbl_dir, "allBact_toFamilyID.tsv"))

# For the sequences in each family, do pairwise alignment and clustering and select
familySubsets_l <- lapply(1:length(bact_famNames), function(index) {

	famName	<- bact_famNames[index]
	message(famName)

	famSubset	<- getFamSubset(famName)
	return(famSubset)
})
familySubsets_df	<- bind_rows(familySubsets_l)

# Extract and write out just the bacterial sequences in the by-family subset
bact_AAsubset		<- bact_AAset[which(names(bact_AAset) %in% familySubsets_df$ID)]
# Reorder sequences by their order in the dataframe
bact_AAsubset		<- bact_AAsubset[order(match(names(bact_AAsubset), familySubsets_df$ID))]

writeXStringSet(bact_AAsubset, filepath = file.path(bact_dir, "HU_IHF_bacteria_famSubset.fasta"))
write_tsv(familySubsets_df, path = file.path(tbl_dir, "byFamSubsetBact_select.tsv"))


## Process the domains for the subset of bacteria
fasta	<- read.fasta(file = file.path(bact_dir, "HU_IHF_bacteria_famSubset.fasta"), seqtype = "AA")
doms	<- read_tsv(file = file.path(dom_dir, "domains_bacteriaSubset.txt"), comment = "")

# Process the domains data
# wDoms	<- processDomains(doms, familySubsets_df)

domains_fix <- doms %>% 
	mutate(Query = str_trim(str_extract(Query, ">+.+?$"))) %>%
	mutate(Query = str_replace(Query, ">", "")) %>%
	dplyr::rename(ID = Query) %>%
	unite(combPos, From, To, Incomplete, sep = " ", remove = FALSE) %>%
	mutate(combPos = paste0("[", combPos, "]"))

by_nDomain	<- domains_fix %>%
	group_by(ID) %>%
	summarise(numDoms = n(), domPos = paste(combPos, collapse = ";"))
combPos	<- left_join(familySubsets_df, by_nDomain, by = "ID")

absentDomIDs	<- combPos %>% filter(is.na(numDoms)) %>% pull(ID)
combPos			<- combPos %>% filter(!is.na(numDoms))
fasta			<- fasta[which(!names(fasta) %in% absentDomIDs)]

# For the split data frame
splitPos	<- domains_fix %>%
	select(ID, combPos, From, To, Incomplete) %>%
	left_join(familySubsets_df, by = "ID") %>%
	group_by(ID) %>%
	arrange(From) %>%
	mutate(Section = case_when(
		n() == 1 ~ "Full",
		From == min(From) ~ "Nterm",
		TRUE ~ "Cterm")) %>%
	mutate(newID = paste(newName, Section, sep = "_")) %>%
	mutate(newID = str_replace_all(newID, " ", "_")) %>%
	arrange(ID)


# Write out the reformatted fasta - not split by domain (but containing domain info)
header <- combPos %>% 
	unite(Heading, newName, ID, Function, Species, taxid, numDoms, domPos, sep = " | ", remove = FALSE)
new_header	<- setNewHeader(fasta, header)
write.fasta(new_header, names = names(new_header), file.out = file.path(bact_dir, paste0("HU_IHF_bacteriaSubset_wTaxid.fasta")))

# Split the sequences containing multiple domains
splitByDom	<- splitFastaByDomain(fasta, splitPos)
splitHeader <- splitPos %>% 
	unite(Heading, newID, ID, Function, Species, taxid, combPos, sep = " | ", remove = FALSE)

## Set new names for the fasta
splitNewHd <- mclapply(1:length(splitByDom), function(seqIndex) {

	seqName		<- names(splitByDom)[seqIndex]
	seq			<- splitByDom[[seqIndex]]
	attrName	<- attributes(seq)$name

	if (!identical(seqName, attrName)) stop("ERROR: object and attribute names are not the same")

	# Get the relevant (updated) header, and assign it as the "Annot" attribute to the sequence
	matchHeader <- splitHeader %>% filter(newID == seqName) %>% pull(Heading)
	attributes(seq)$Annot	<- matchHeader

	return(seq)
}, mc.cores = 10)

newNames	<- lapply(splitNewHd, function(seq) {
	name	<- attributes(seq)$Annot
	return(name)
})
names(splitNewHd)	<- newNames

## ------------------ ##
# Write out the various sequences to fasta

# We will split into single domains (complete), single domains (incomplete), and double domains (all)
singleDoms_IDs	<- combPos %>% filter(numDoms == 1) %>% pull(ID)
singleDom_cmpl	<- splitHeader %>% filter(ID %in% singleDoms_IDs & Incomplete == "-") %>% pull(Heading)
singleDom_incm	<- splitHeader %>% filter(ID %in% singleDoms_IDs & Incomplete != "-") %>% pull(Heading)
multiDom		<- splitHeader %>% filter(!ID %in% singleDoms_IDs) %>% pull(Heading)

# Write the split fasta files
singleDom_cmpl_fasta	<- splitNewHd[which(names(splitNewHd) %in% singleDom_cmpl)]
singleDom_incm_fasta	<- splitNewHd[which(names(splitNewHd) %in% singleDom_incm)]
multiDom_fasta			<- splitNewHd[which(names(splitNewHd) %in% multiDom)]

write.fasta(splitNewHd, names = names(splitNewHd), file.out = file.path(bact_dir, "HU_IHF_bacteriaSubset_byDomSplit_all.fasta"))
write.fasta(singleDom_cmpl_fasta, names = names(singleDom_cmpl_fasta), file.out = file.path(bact_dir, "HU_IHF_bacteriaSubset_oneDom_complete.fasta"))
write.fasta(singleDom_incm_fasta, names = names(singleDom_incm_fasta), file.out = file.path(bact_dir, "HU_IHF_bacteriaSubset_oneDom_incomplete.fasta"))
write.fasta(multiDom_fasta, names = names(multiDom_fasta), file.out = file.path(bact_dir, "HU_IHF_bacteriaSubset_twoDom_all.fasta"))

## ------------------ ##
# Write out the various sequences to fasta with reduced headers for alignment

forAlignBact_dir	<- file.path(bact_dir, "fastaForAlign")
if (!dir.exists(forAlignBact_dir)) dir.create(forAlignBact_dir)

write.fasta(splitNewHd, names = str_trim(str_split_fixed(names(splitNewHd), "\\|", 2)[,1]), file.out = file.path(forAlignBact_dir, "HU_IHF_bacteriaSubset_byDomSplit_all.fasta"))
write.fasta(singleDom_cmpl_fasta, names = str_trim(str_split_fixed(names(singleDom_cmpl_fasta), "\\|", 2)[,1]), file.out = file.path(forAlignBact_dir, "HU_IHF_bacteriaSubset_oneDom_complete.fasta"))
write.fasta(singleDom_incm_fasta, names = str_trim(str_split_fixed(names(singleDom_incm_fasta), "\\|", 2)[,1]), file.out = file.path(forAlignBact_dir, "HU_IHF_bacteriaSubset_oneDom_incomplete.fasta"))
write.fasta(multiDom_fasta, names = str_trim(str_split_fixed(names(multiDom_fasta), "\\|", 2)[,1]), file.out = file.path(forAlignBact_dir, "HU_IHF_bacteriaSubset_twoDom_all.fasta"))

## ------------------ ##
# Separate very long sequences

shortComplete	<- singleDom_cmpl_fasta[which(lapply(singleDom_cmpl_fasta, length) <= 110)]
longComplete	<- singleDom_cmpl_fasta[which(lapply(singleDom_cmpl_fasta, length) > 110)]

write.fasta(shortComplete, names = str_trim(str_split_fixed(names(shortComplete), "\\|", 2)[,1]), file.out = file.path(forAlignBact_dir, "HU_IHF_bacteriaSubset_oneDom_complete_short.fasta"))
write.fasta(longComplete, names = str_trim(str_split_fixed(names(longComplete), "\\|", 2)[,1]), file.out = file.path(forAlignBact_dir, "HU_IHF_bacteriaSubset_oneDom_complete_long.fasta"))



####

# mafft-linsi --thread 10 --localpair HU_IHF_bacteriaSubset_oneDom_complete_short.fasta > ../align/shortLINSI.fasta
# mafft-linsi --thread 10 --localpair --add ../../Archaea/HU_IHF_archaea_oneDom_complete.fasta  --reorder ../align/shortLINSI.fasta > ../align/shortAddArch.fasta
# t_coffee HU_IHF_bacteriaSubset_oneDom_complete_short.fasta -mode expresso -output=fasta_aln -outfile=../align/shortExpresso.fasta
# mafft-linsi --thread 10 --localpair --add ../../Archaea/HU_IHF_archaea_oneDom_complete.fasta  --reorder ../align/shortExpresso.fasta > ../align/shortExpressoAddArch.fasta

# mafft-linsi --thread 10 --localpair --add  ../../../Archaea/HU_IHF_archaea_oneDom_complete.fasta  --reorder bactShort_cobaltOnlineAln.fa > bactShort_cobaltOnlineAln_addArch.fa

exec raxml -f a -s /Users/aesin/Desktop/HTa_tonio/Round2/Align/Combined/COBALT/comb_allOneDom_complete_aln_relab.fa -n complete_aln_relab_raxml_tree.txt -m PROTCATAUTO -p $RANDOM -x $RANDOM -N 100 -T 20




