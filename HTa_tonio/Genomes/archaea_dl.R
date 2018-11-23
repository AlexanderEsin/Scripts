library(tidyverse)

## Set up directories
genome_dir		<- "/users/aesin/Desktop/HTa_tonio/Genomes/archaea"
genome_out_dir	<- file.path(genome_dir, "Genome_gbffs")
genome_list_dir	<- file.path(genome_dir, "Genome_lists")


# Create dirs if they don't exist
if (!dir.exists(genome_out_dir)) dir.create(genome_out_dir)
if (!dir.exists(genome_list_dir)) dir.create(genome_list_dir)


if (!file.exists(file.path(genome_list_dir, "All_arch_genomes.tsv"))) {
	# Open the RefSeq annotation summary table
	genome_tbl <- read.table(file.path(genome_dir, "assembly_summary_151018.txt"), sep = "\t", as.is = TRUE, header = TRUE, fill = TRUE, quote = "", stringsAsFactors = FALSE, comment.char = "")

	# Find unique taxids
	unique_taxids			<- unique(genome_tbl$taxid)

	# Filter out all genomes with a single completely assembled genomes
	count_assemblies		<- as.data.frame(table(genome_tbl$taxid), stringsAsFactors = FALSE)
	single_ass_taxids		<- as.numeric(count_assemblies$Var1[count_assemblies$Freq == 1])
	single_ass_df			<- genome_tbl[which(genome_tbl$taxid %in% single_ass_taxids),]

	# Taxids represented by more than one assembly
	multi_complete_ass		<- genome_tbl[!(genome_tbl$taxid %in% single_ass_taxids),]
	multi_ass_taxids		<- unique(multi_complete_ass$taxid)

	# Assembly completeness levels
	completeness_lvl		<- data.frame(assembly_level = c("Complete Genome", "Chromosome", "Scaffold", "Contig"), stringsAsFactors = FALSE)

	multi_select_list	<- lapply(multi_ass_taxids, function(taxid) {
		subset_data	<- multi_complete_ass[which(multi_complete_ass$taxid == taxid),]

		# Check if one is representative
		if (any(subset_data$refseq_category == "representative genome")) {
			select_ass	<- subset_data[which(subset_data$refseq_category == "representative genome"),]
			if (nrow(select_ass) == 1) return(select_ass) else {
				message(paste0("Multiple representative genomes at taxid ", taxid, " - returning first entry"))
				return(select_ass[1,])
			}
		}

		# Pick the best assembled
		order_by_ass_lvl	<- left_join(completeness_lvl, subset_data, by = "assembly_level")
		# Remove any completeness levels without any genomes
		order_by_ass_lvl	<- order_by_ass_lvl[!is.na(order_by_ass_lvl$assembly_accession),]
		# Select the top represented assembly level
		highest_assembly	<- order_by_ass_lvl[1,1]
		# Select all assemblies at the highest assembly level
		all_highest	<- order_by_ass_lvl %>% filter(assembly_level == highest_assembly)

		# If there is only one assembly at this level, return it
		if (nrow(all_highest) == 1) return(all_highest)

		# Otherwise sort by latest
		latest_date <- max(as.Date(all_highest$seq_rel_date))
		# Pick the assembly with the latest date...
		latest_assembly	<- all_highest %>% filter(seq_rel_date == latest_date)
		# If there's only one with the latest date, return it
		if (nrow(latest_assembly) == 1) return(latest_assembly)

		# Otherwise pick one randomly...
		random_index	<- sample(1:nrow(latest_assembly), 1)
		random_assembly	<- latest_assembly[random_index,]

		return(random_assembly)
	})

	# Bind the selected entries which had multiple assemblies into df
	multi_select_df	<- bind_rows(multi_select_list)

	# Bind the single assembly and multi-assembly dataframes together
	all_select_df	<- bind_rows(multi_select_df, single_ass_df)

	# Write out the full table of downloaded genomes
	write.table(all_select_df, file = file.path(genome_list_dir, "All_arch_genomes.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

	# Write out just a list of the acc_ass_names
	acc_ass_names	<- lapply(all_select_df$assembly_accession, function(accession) {
		name	<- paste(c(accession, all_select_df$asm_name[which(all_select_df$assembly_accession == accession)]), collapse = "_")
		return(name)
	})
	write.table(unlist(acc_ass_names), file = file.path(genome_list_dir, "Acc_ass_list.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)

} else {
	all_select_df	<- read_tsv(file.path(genome_list_dir, "All_arch_genomes.tsv"), col_names = TRUE, quote = "")
}



# Foreach completely assembled genome chosen above, 
# download the gbff file with aria2c.
acc_ass_tax <- lapply(1:nrow(all_select_df), function(row_index) {
	genome_info	<- all_select_df[row_index,]
	genome_name	<- genome_info$organism_name
	ftp_path	<- genome_info$ftp_path
	taxid		<- genome_info$taxid

	## Acc_ass will be used to ID the files
	acc_ass		<- basename(ftp_path)
	file_name	<- paste0(acc_ass, "_genomic.gbff.gz")
	dl_file		<- file.path(ftp_path, file_name)

	if (!file.exists(file.path(genome_out_dir, file_name))) {
		message(paste0("Dowloading... ", row_index, " // ", nrow(all_select_df)))
		system2("aria2c", args = c(paste0("--dir ", genome_out_dir), paste0("--out ", file_name), "-q", dl_file))
	}
	return(data.frame(Acc_ass = acc_ass, Taxid = taxid, stringsAsFactors = FALSE))
})

# Write out a table of taxids and acc_ass IDs
acc_ass_dled_df	<- bind_rows(acc_ass_tax)
write.table(acc_ass_dled_df, file = file.path(genome_list_dir, "Acc_ass_taxid_table.tsv"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)



## There are 5073 downloaded genomic gbff files
dled_genome_num	<- length(dir(genome_out_dir))