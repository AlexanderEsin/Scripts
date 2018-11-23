## 12/11/2017
## Filter out all the "complete genomes" and "chromosomes" from the Bacterial Refseq db: ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/
## Downloaded assembly_summary_refseq.txt from above url: wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt; mv assembly_summary.txt assembly_summary_refseq_131117.txt
## Manually remove the hashtag infront of columns names (line 2)

library(tidyverse)

## Set up directories
genome_dir		<- "/users/aesin/Desktop/Geo_again/Genomes"
genome_out_dir	<- file.path(genome_dir, "Genome_gbffs")
genome_list_dir	<- file.path(genome_dir, "Genome_lists")

# Create dirs if they don't exist
if (!dir.exists(genome_out_dir)) dir.create(genome_out_dir)
if (!dir.exists(genome_list_dir)) dir.create(genome_list_dir)

if (!file.exists(file.path(genome_list_dir, "All_complete_genomes.tsv"))) {
	
	## Open the RefSeq annotation summary table
	bact_genome_tbl <- read.table(file.path(genome_dir, "assembly_summary_refseq_131117.txt"), sep = "\t", as.is = TRUE, header = TRUE, fill = TRUE, quote = "", stringsAsFactors = FALSE, comment.char = "")

	## Select all genomes which are 'complete' - this is 8296 genomes
	bact_genome_complete	<- bact_genome_tbl[which(bact_genome_tbl$assembly_level == "Complete Genome"),]
	## These 8296 genomes represent 5073 unique taxonomic IDs (i.e. different species)
	unique_taxids			<- unique(bact_genome_complete$taxid)

	## Filter out all genomes with a single completely assembled genomes
	count_assemblies		<- as.data.frame(table(bact_genome_complete$taxid), stringsAsFactors = FALSE)
	single_ass_taxids		<- as.numeric(count_assemblies$Var1[count_assemblies$Freq == 1])
	single_complete_ass_df	<- bact_genome_complete[which(bact_genome_complete$taxid %in% single_ass_taxids),]

	## 4577 / 5073 taxids have a single assembly only
	## For the remainder of taxids (496) we select the latest complete assembly.
	## If > 1 assembly has the same (latest) date, select one at random
	multi_complete_ass		<- bact_genome_complete[!(bact_genome_complete$taxid %in% single_ass_taxids),]
	multi_ass_taxids		<- unique(multi_complete_ass$taxid)

	selected_multi_ass	<- lapply(multi_ass_taxids, function(taxid) {
		## Get all the entries corresponding to the taxid (must be > 1)
		subset_data	<- multi_complete_ass[which(multi_complete_ass$taxid == taxid),]
		latest_date <- max(as.Date(subset_data$seq_rel_date))
		
		## Pick the assembly with the latest date
		picked_assembly	<- subset_data[which(subset_data$seq_rel_date == latest_date),]

		paste0("\n")

		## If there's more than one, pick one at random
		if (nrow(picked_assembly) > 1) {
			random_index	<- sample(1:nrow(picked_assembly), 1)
			picked_assembly	<- picked_assembly[random_index,]
			picked_randomly <- TRUE
		} else {
			picked_randomly <- FALSE
		}

		## Return a dataframe (all columns as in bact_genome_tbl) + a new column
		## indicating which assemblies had to be selected randomly
		selected_assembly	<- cbind(picked_assembly, random_pick = picked_randomly)
		return(selected_assembly)
	})

	selected_multi_ass_df	<- bind_rows(selected_multi_ass)

	## 136 taxids had the assembly chosen randomly 
	single_complete_ass_df$random_pick <- FALSE
	all_complete_ass_bact	<- rbind(single_complete_ass_df, selected_multi_ass_df)

	## 25 Anoxybacillus and (para)Geobacillus genomes in this dataset.
	## Write these out into a seperate table
	matchAG			<- paste(c("Geobacillus", "Anoxybacillus"), collapse = "|")
	AG_in_dataset	<- all_complete_ass_bact[grep(matchAG, all_complete_ass_bact$organism_name, ignore.case = TRUE),]


	## Write out the lists of genomes (one containing all the
	## genomes downloaded, the other only the AG genomes). Also
	genome_list_dir	<- file.path(genome_dir, "Genome_lists")
	if (!dir.exists(genome_list_dir)) dir.create(genome_list_dir)
	# write.table(all_complete_ass_bact, file = file.path(genome_list_dir, "All_complete_genomes.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
	# write.table(AG_in_dataset, file = file.path(genome_list_dir, "AG_complete_genomes.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

	## Write out a column of just acc_ass for the AG genomes.
	acc_ass_names	<- lapply(AG_in_dataset$assembly_accession, function(accession) {
		name	<- paste(c(accession, AG_in_dataset$asm_name[which(AG_in_dataset$assembly_accession == accession)]), collapse = "_")
		return(name)
	})
	# write.table(unlist(acc_ass_names), file = file.path(genome_list_dir, "AG_acc_ass_names.txt"), sep = "\n", row.names = FALSE, col.names = FALSE, quote = FALSE)

} else {
	# all_complete_ass_bact	<- read.table(file.path(genome_list_dir, "All_complete_genomes.tsv"), sep = "\t", as.is = TRUE, header = TRUE, fill = TRUE, quote = "", stringsAsFactors = FALSE, comment.char = "")
}


## Foreach completely assembled genome chosen above, 
## download the gbff file with aria2c.
acc_ass_tax <- lapply(1:nrow(all_complete_ass_bact), function(row_index) {
	genome_info	<- all_complete_ass_bact[row_index,]
	genome_name	<- genome_info$organism_name
	ftp_path	<- genome_info$ftp_path
	taxid		<- genome_info$taxid

	## Acc_ass will be used to ID the files
	acc_ass		<- basename(ftp_path)
	file_name	<- paste0(acc_ass, "_genomic.gbff.gz")
	dl_file		<- file.path(ftp_path, file_name)

	if (!file.exists(file.path(genome_out_dir, file_name))) {
		message(paste0("Dowloading... ", row_index, " // ", nrow(all_complete_ass_bact)))
		system2("aria2c", args = c(paste0("--dir ", genome_out_dir), paste0("--out ", file_name), "-q", dl_file))
	}
	return(data.frame(Acc_ass = acc_ass, Taxid = taxid, stringsAsFactors = FALSE))
})

## Write out a table of taxids and acc_ass IDs
acc_ass_dled_df	<- bind_rows(acc_ass_tax)
# write.table(acc_ass_dled_df, file = file.path(genome_list_dir, "Acc_ass_taxid_table.tsv"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)



## There are 5073 downloaded genomic gbff files
dled_genome_num	<- length(dir(genome_out_dir))







## ---- // ---- ##
## Snippet below used to recover the final genome tables
## they were written with the wrong separators initially

# gbff_dled_list	<- Sys.glob(file.path(genome_out_dir, "*.gz"))
# gbff_base		<- basename(gbff_dled_list)
# gbff_clean		<- str_replace(gbff_base, "_genomic.gbff.gz", "")
# complete_genomes_list	<- lapply(gbff_clean, function(gbff) {
# 	print(gbff)
# 	found_row	<- bact_genome_tbl[grep(gbff, bact_genome_tbl$ftp_path, ignore.case = TRUE),]
# 	if (nrow(found_row) == 1) {
# 		return(found_row)
# 	} else {
# 		warning("Either found 0 or > 1 rows")
# 	}
# })
# complete_genomes_df	<- bind_rows(complete_genomes_list)
# matchAG			<- paste(c("Geobacillus", "Anoxybacillus"), collapse = "|")
# AG_in_dataset	<- complete_genomes_df[grep(matchAG, complete_genomes_df$organism_name, ignore.case = TRUE),]



# write.table(AG_in_dataset, file = file.path(genome_list_dir, "AG_complete_genomes2.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
# write.table(complete_genomes_df, file = file.path(genome_list_dir, "All_complete_genomes2.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
