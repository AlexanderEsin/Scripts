#!/usr/bin/tclsh

## In this script we create an sqlite database to hold
## the bitscores for all orthologs, we will comparalog
## the bitscores of possible in-paralogs to ortholog scores

source /home/ade110/Scripts/General_utils.tcl
package require sqlite3

set evalue		[lindex $argv 0]
set eval_trunc	[string range $evalue 2 end]


## Set relevant directories and sqlite DB name
set direct			/scratch/ade110/Geo_again/RBBH

set para_comp_dir	$direct/Paralog_comparison
set comparison_db	$direct/RBBH_OrthCompareDB
set db_build_log	$direct/RBBH_OrthCompareDB.log

## If the paralog comparison DB does not exist, create
if {[file exists $comparison_db] != 1} {
	## Open log
	set logchan			[open $db_build_log w]
	puts $logchan		"Beginning to build database from RBBHs at evalue: $eval_trunc"

	## Make an sqliteDB of the ortholog RBBHs with bitscores.
	## We will use this as a comparison set
	sqlite3 db1 $comparison_db
	db1 eval {PRAGMA main.page_size = 4096}
	db1 eval {PRAGMA main.cache_size=10000}
	db1 eval {PRAGMA main.locking_mode=EXCLUSIVE}
	db1 eval {PRAGMA main.synchronous=NORMAL}
	db1 eval {PRAGMA main.journal_mode=WAL}

	db1 eval {
		create table t1(GeneA text, GeneB text, bitscore integer)
	}

	## List of all the RBBH scores (bitscore) against which
	## we will try to indentify in-paralogs
	set RBBH_comp_files [glob $para_comp_dir/*txt]
	set num_comp_files	[llength $RBBH_comp_files]
	set counter 1

	foreach RBBH_comp_file $RBBH_comp_files {
		## Data is a list of all the bitscore RBBH comparisons
		## Each line is formatted: GeneA GeneB Bitscore
		set data		[split [string trim [openfile $RBBH_comp_file]] \n]

		## Do not include the header
		set no_header 	[lrange $data 1 end]

		foreach RBBH $data {
			set values	[split $RBBH \t]
			set geneA	[lindex $values 0]
			set geneB	[lindex $values 1]
			set bits	[lindex $values 2]

			db1 eval {
				insert into t1 values($geneA,$geneB,$bits)
			}
		}

		puts $logchan "Adding files to DB: $counter / $num_comp_files ..."
		incr counter
	}
	puts $logchan "\nAdding files to DB: DONE"
	puts $logchan "Making database indexes ..."
	## Create indexes on both the gene columns
	db1 eval {CREATE INDEX old_id_index ON t1 (GeneA)}
	db1 eval {CREATE INDEX new_id_index ON t1 (GeneB)}
	puts $logchan "Making database indexes: DONE"
	close $logchan

## Otherwise just open the database
} else {
	puts stdout "Database already exists!"
	sqlite3 db1 $comparison_db
}
