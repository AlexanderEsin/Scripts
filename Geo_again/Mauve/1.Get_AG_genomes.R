master_path		<- "/Users/aesin/Desktop/Geo_again"
genome_path		<- file.path(master_path, "Genomes")
genomeL_path	<- file.path(genome_path, "Genome_lists")
gbffs_path		<- file.path(genome_path, "Genome_gbffs")

AG_gbff_path	<- file.path(genome_path, "AG_genome_gbffs")
if (!dir.exists(AG_gbff_path)) dir.create(AG_gbff_path)

AG_accAss_fl	<- file.path(genomeL_path, "AG_acc_ass_names.txt")

AG_accAss_tbl	<- read.table(file = AG_accAss_fl, sep = "\n", header = FALSE, stringsAsFactors = FALSE)
AG_accAss_l		<- AG_accAss_tbl$V1

AG_gbff_files	<- file.path(gbffs_path, paste0(AG_accAss_l, "_genomic.gbff.gz"))

# Check files exist
filesPresent	<- ifelse(all(unlist(lapply(AG_gbff_files, file.exists))), TRUE, FALSE)

if (filesPresent) {
	message("All files found")
	invisible(lapply(AG_gbff_files, file.copy, to = AG_gbff_path))
	message("All files copied")
} else {
	stop("Not all files found in source directory")
}