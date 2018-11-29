#!/usr/bin/env Rscript

# Any libraries explicitly used in the script
require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load(tidyverse, RSQLite, fs)

master			<- file.path("/Users/aesin/Desktop/Staph")
grp_list_dir	<- file.path(master, "Staph_genomes", "Genome_lists")
all_db_path		<- file.path(master, "All_prot_db_new")
grp_db_path		<- file.path(master, "ForStaph_prot_db")

toRemove		<- read_tsv(file.path(grp_list_dir, "grpToRemove.tsv"))

# Copy the database
if (!file.exists(grp_db_path)) file_copy(all_db_path, grp_db_path)


# ----------------------------------------------------- #
## Drop the is_ag, COGcat, and OrthGroup columns
## // IN SHELL WITHIN SQLITE RUN:

# PRAGMA main.page_size = 4096;
# PRAGMA main.cache_size=10000;
# PRAGMA main.locking_mode=EXCLUSIVE;
# PRAGMA main.synchronous=NORMAL;
# PRAGMA main.journal_mode=WAL;

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
# VACUUM;
# .exit

# ----------------------------------------------------- #
# Open database
dbConn			<- dbConnect(RSQLite::SQLite(), grp_db_path)

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
