#!/usr/bin/env Rscript

# Any libraries explicitly used in the script
require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load(tidyverse, RSQLite, fs)

master			<- file.path("/Users/aesin/Desktop/Bacillus")
bac_list_dir	<- file.path(master, "Bac_genomes", "Genome_lists")
all_db_path		<- file.path(master, "All_prot_db_new")
bac_db_path		<- file.path(master, "ForBac_prot_db")

toRemove		<- read_tsv(file.path(bac_list_dir, "bacillusToRemove.tsv"))

# Copy the database
if (!file.exists(bac_db_path)) file_copy(all_db_path, bac_db_path)


# ----------------------------------------------------- #
## Drop the is_ag, COGcat, and OrthGroup columns
## // IN SHELL WITHIN SQLITE RUN:

# PRAGMA table_info(t1);

# CREATE TABLE new_t1 (
# 	protID text,
# 	locus text,
# 	gene_start integer,
# 	gene_end integer,
# 	strand text,
# 	plasmid text,
# 	product text,
# 	sequence text,
# 	NuclSeq text,
# 	acc_ass text,
# 	binomial text,
# 	strain text,
# 	taxid integer,
# 	genome_l integer
# );
 
# INSERT INTO new_t1 SELECT protID, locus, gene_start, gene_end, strand, plasmid, product, sequence, NuclSeq, acc_ass, binomial, strain, taxid, genome_l FROM t1;

# DROP TABLE IF EXISTS t1; 
# ALTER TABLE new_t1 RENAME TO t1;

# CREATE INDEX protID_index on t1 (protID);
# CREATE INDEX taxid_index on t1 (taxid);
# CREATE INDEX locTag_index on t1 (locus);
# .exit

# ----------------------------------------------------- #
# Open database
dbConn			<- dbConnect(RSQLite::SQLite(), bac_db_path)

## Remove all rows corresponding to the bacillus species that are no longer wanted
dbBegin(dbConn)
rs <- dbSendStatement(dbConn, "DELETE FROM t1 WHERE taxid = :taxids")
dbBind(rs, param = list(taxids = toRemove$Taxid))
dbClearResult(rs)
dbCommit(dbConn)

## Confirm the right number of taxids were removed
# query	<- dbSendQuery(dbConn, 'SELECT COUNT (DISTINCT taxid) FROM t1')
# count	<- dbFetch(query)
# dbClearResult(query)

dbDisconnect(dbConn)
