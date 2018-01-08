#!/usr/local/bin/tclsh

source ~/Documents/Scripts/General_utils.tcl
package require sqlite3

set group_name	[lindex $argv 0]
set evalue		[lindex $argv 1]
set eval_trunc	[string range $evalue 2 end]


## Set relevant directories and sqlite DB name
set direct			/users/aesin/Desktop/Geo_again/Consensus_groups/$group_name
set para_comp_dir	$direct/RBBH/Paralog_comparison
set para_poss_dir	$direct/Find_InParalogs
set para_RBBH_dir	$direct/RBBH/InParalogs_RBBH

file mkdir			$para_RBBH_dir

set comparison_db	$direct/RBBH/RBBH_paralog_comparison_db


## If the paralog comparison DB does not exist, create
if {[file exists $comparison_db] != 1} {
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

		puts -nonewline stdout "Adding files to DB: $counter / $num_comp_files ...\r"
		flush stdout
		incr counter
	}
	puts "\nAdding files to DB: DONE"
	puts "Making database indexes ..."
	## Create indexes on both the gene columns
	db1 eval {CREATE INDEX old_id_index ON t1 (GeneA)}
	db1 eval {CREATE INDEX new_id_index ON t1 (GeneB)}
	puts "Making database indexes: DONE"

## Otherwise just open the database
} else {
	puts "Database already exists!"
	sqlite3 db1 $comparison_db
}


## Sets of possible paralogs are derived from:
## Scripts/Geo_again/Anogeo/Find_paralogs/Find_paralogs.tcl
set poss_paralog_files	[glob $para_poss_dir/Poss_paralogs/*tsv]

set true_paralog_list	{"Paralog_A\tParalog_B"}
set paralog_for_RBBH	{}
set paralog_counter		0
set process_counter		1

puts "Identifying true paralogs ..."

## Possible paralogs are split by species
foreach paralog_file $poss_paralog_files {

	## Read in the possible paralogs data
	set poss_para_data	[split [string trim [openfile $paralog_file]] \n]

	puts -nonewline stdout "Processing paralog file: $process_counter // [llength $poss_paralog_files]\r"
	flush stdout
	incr process_counter

	foreach possible_paralogs $poss_para_data {

		## Iterate over each possible paralogous combination
		set values		[split $possible_paralogs \t]
		set para_A		[lindex $values 0]
		set para_B		[lindex $values 1]
		set para_bits	[lindex $values 2]
		set para_eval	[lindex $values 3]

		## Find the maximum ortholog bitscore for either potential paralog
		## against any other orthologous gene
		set orth_bits [db1 eval {
			select max(bitscore) from t1 where GeneA = $para_A or GeneB = $para_A or GeneA = $para_B or GeneB = $para_B
		}]

		## If there is no value, the above command returns "{}"
		## We treat that as a list with an empty first element
		set orth_bits [lindex $orth_bits 0]
		if {[string length $orth_bits] == 0} {
			set orth_bits 0
		}

		## If the paralog bitscore is the same or higher than the
		## top bitscore against any ortholog and the evalue between
		## the paralogs is below the treshold, we add this paralogy
		## to a list of true paralogs.
		if {$para_bits >= $orth_bits && $para_eval <= $evalue} {
			lappend true_paralog_list	"$para_A\t$para_B"

			lappend paralog_for_RBBH	"$para_A\t$para_B\t$para_eval"
			lappend paralog_for_RBBH	"$para_B\t$para_A\t$para_eval"
			incr paralog_counter	
		}

	}
}
puts "\nIdentifying true paralogs: DONE"
puts "Identified $paralog_counter true paralogs"

## Write out the list of true paralogs
set out		[open $para_poss_dir/True_paralogs\_$eval_trunc.tsv w]
puts $out	[join $true_paralog_list \n]
close $out

puts "Adding paralog RBBHs to RBBH directory ..."

## Write out the InParalog specific RBBH file into the
## RBBH directory
set out		[open $para_RBBH_dir/InParalogs_RBBH_$eval_trunc\.txt w]
puts $out	[join $paralog_for_RBBH \n]
close $out

puts "Added [expr [llength $paralog_for_RBBH] / 2] InParalog RBBHs at $evalue to:\n\t$direct/RBBH/InParalogs_RBBH_$eval_trunc\.txt"

db1 close
