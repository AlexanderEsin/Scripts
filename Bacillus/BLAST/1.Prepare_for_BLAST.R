# Any libraries explicitly used in the script
require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load(tidyverse, RSQLite, fs)

# ----------------------------------------------------- #
master			<- file.path("/Users/aesin/Desktop/Bacillus")
blast_dir		<- file.path(master, "BLAST")

# BLAST sub directories
blast_inAll		<- file.path(blast_dir, "IN_all")
blast_inSplit	<- file.path(blast_dir, "IN_all_split")
blast_dbAll		<- file.path(blast_dir, "DB_all")
blast_inGrp		<- file.path(blast_dir, "IN_grp")
blast_dbGrp		<- file.path(blast_dir, "DB_grp")

# Create the directories
dir_create(c(blast_inAll, blast_inSplit, blast_dbAll, blast_inGrp, blast_dbGrp))

# Proteome directory
clean_prot_dir	<- file.path(master, "Proteomes", "Clean")

# List of the group accessions (against which we BLAST)
bac_list_dir	<- file.path(master, "Bac_genomes", "Genome_lists")
toKeep_tbl		<- read_tsv(file.path(bac_list_dir, "bacillusToKeep.tsv"))
accAss_tax_tbl	<- read_tsv(file.path(bac_list_dir, "bac_AccAssTaxid_table.tsv"))

# ----------------------------------------------------- #
invisible(lapply(1:nrow(accAss_tax_tbl), function(index) {

	acc_ass		<- accAss_tax_tbl[index,] %>% pull(Acc_ass)
	fasta_file	<- file.path(clean_prot_dir, paste0(acc_ass, "_proteome.fasta.gz"))
	db_name		<- file.path(blast_dbAll, acc_ass)

	# Copy file into the IN_all directory
	file_copy(fasta_file, blast_inAll)

	# Make the database
	system2("makeblastdb", args = c("-dbtype prot", paste0("-out ", db_name), paste0("-title ", acc_ass)), input = read_lines(gzfile(fasta_file)), stdout = NULL)

	# If this accession is one of the core group, copy accordingly
	if (acc_ass %in% toKeep_tbl$Acc_ass) {
		# Copy into IN_grp
		file_copy(fasta_file, blast_inGrp)
		# Make db
		grp_db_name	<- file.path(blast_dbGrp, acc_ass)
		system2("makeblastdb", args = c("-dbtype prot", paste0("-out ", grp_db_name), paste0("-title ", acc_ass)), input = read_lines(gzfile(fasta_file)), stdout = NULL)
		# Ticker
		message(paste0("--> Group", index))
		return(NULL)
	}
	message(index)
	return(NULL)
}))


# ----------------------------------------------------- #
all_proteomes	<- dir(blast_inAll, full.names = TRUE)
per_folder		<- 3
i = 0

for (f in all_proteomes) {
	subdir	<- file.path(blast_inSplit, paste0("dir_", sprintf("%04d", floor((i / per_folder) + 1))))
	dir_create(subdir)
	file_move(f, subdir)
	i = i + 1
}

dir_delete(blast_inAll)

# ----------------------------------------------------- #




