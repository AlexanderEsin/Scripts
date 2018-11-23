#!/usr/bin/env Rscript

require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load("RSQLite", "tidyverse", "parallel", "seqinr")

master_dir		<- "/Users/aesin/Desktop/HTa_tonio"
blast_out_dir	<- file.path(master_dir, "Blast", "Blast_result")
alignment_dir	<- file.path(master_dir, "Alignment")
# Write out the subsets of sequences
if (!dir.exists(alignment_dir)) dir.create(alignment_dir)


# Open All_prot database
HTa_db_path		<- file.path(master_dir, "HTa_prot_db")
dbConn			<- dbConnect(RSQLite::SQLite(), HTa_db_path)



# Get a list of the blast output files, and see how many have hits
all_out_files	<- dir(blast_out_dir, full.names = TRUE)
total_blast_num	<- length(all_out_files)
blast_w_hits	<- all_out_files[which(lapply(all_out_files, file.size) > 0)]

# Extract the blast results from each file with a hit
blastOut_score_list	<- mclapply(blast_w_hits, function(file) {
	data <- read_tsv(file, col_names = c("query", "target", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send",
   "evalue", "bitscore"), col_types = cols())
	return(data)
}, mc.cores = 10)

# Combine into data frame, remove the query column
blastOut_score_df	<- bind_rows(blastOut_score_list) %>% select(-query) %>% arrange(evalue)

# bitscore distribution
ggplot(data = blastOut_score_df %>% arrange(desc(bitscore)), aes(x = seq(1:nrow(blastOut_score_df)), y = bitscore)) +
	scale_y_continuous(limits = c(40, 180)) +
	geom_line() +
	theme_classic()

# % ident distribution
ggplot(data = blastOut_score_df %>% arrange(desc(pident)), aes(x = seq(1:nrow(blastOut_score_df)), y = pident)) +
	scale_y_continuous(limits = c(0, 100)) +
	geom_line() +
	theme_classic()

# Histogram of log evalues
ggplot(data = blastOut_score_df, aes(x = log(evalue))) +
	# scale_y_continuous(limits = c(40, 180)) +
	scale_x_continuous(limits = c(-150, 0), breaks = seq(-150, 0, by = 10)) +
	geom_histogram() +
	theme_classic()


# Select a subset for two stage alignemnt
highHomology	<- blastOut_score_df %>% filter(evalue < 1e-15)
lowHomology		<- blastOut_score_df %>% filter(evalue >= 1e-15)


# Extract the fasta sequences

# High homology
protSeqs_tbl	<- dbSendQuery(dbConn, 'SELECT protID, sequence FROM t1 WHERE protID = :protIDs')
dbBind(protSeqs_tbl, param = list(protIDs = highHomology$target))
protSeqs_df	<- dbFetch(protSeqs_tbl)
dbClearResult(protSeqs_tbl)

fastaSeqs <- lapply(protSeqs_df$sequence, as.SeqFastaAA)
write.fasta(sequences = fastaSeqs, names = protSeqs_df$protID, file.out = file.path(alignment_dir, "highHomology.fasta"))

# Low homology
protSeqs_tbl	<- dbSendQuery(dbConn, 'SELECT protID, sequence FROM t1 WHERE protID = :protIDs')
dbBind(protSeqs_tbl, param = list(protIDs = lowHomology$target))
protSeqs_df	<- dbFetch(protSeqs_tbl)
dbClearResult(protSeqs_tbl)

fastaSeqs <- lapply(protSeqs_df$sequence, as.SeqFastaAA)
write.fasta(sequences = fastaSeqs, names = protSeqs_df$protID, file.out = file.path(alignment_dir, "lowHomology.fasta"))


# Write out full blast results table
binomial_tbl	<- dbSendQuery(dbConn, 'SELECT protID, binomial FROM t1 WHERE protID = :protIDs')
dbBind(binomial_tbl, param = list(protIDs = blastOut_score_df$target))
binomial_df	<- dbFetch(binomial_tbl)
dbClearResult(binomial_tbl)


# New column with just first two words of the binomial name
binomial_df %<>% mutate(binomial_trunc = word(binomial, 1, 2, sep = " "))
# Count the number of duplicate organism names, assign them a counter
binomial_df	%<>% group_by(binomial_trunc) %>% mutate(binom_count = seq(1:n()))
# Combine the binomial and counter so each label can be traced back to protID
binomial_df %<>% mutate(new_tip_label = paste(unlist(lapply(str_split(binomial_trunc, " "), paste, collapse = "_")), binom_count, sep = "_"))

# Join to blast output table
fullResults_table 	<- left_join(blastOut_score_df, binomial_df, by = c("target" = "protID"))
write.table(fullResults_table, file = file.path(blast_out_dir, "..", "fullBlastOutput.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)



# Close sqlite3 connection
dbDisconnect(dbConn)