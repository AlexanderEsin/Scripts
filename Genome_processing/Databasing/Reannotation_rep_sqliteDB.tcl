#!/usr/local/bin/tclsh

source ~/Dropbox/Scripts/General_utils.tcl
package require sqlite3

###
set direct /users/aesin/desktop/Clean_Proteomes
set db_name Reannotation_db
set file_to_sqlite_DB Reannotation_report.tsv
###

cd $direct

## Create the sqlite database within the Rbhh folder ##
sqlite3 db1 $db_name
db1 eval {PRAGMA main.page_size = 4096}
db1 eval {PRAGMA main.cache_size=10000}
db1 eval {PRAGMA main.locking_mode=EXCLUSIVE}
db1 eval {PRAGMA main.synchronous=NORMAL}
db1 eval {PRAGMA main.journal_mode=WAL}
db1 eval {create table t1(old_id text, new_id text, accession text, old_start text, old_end text, old_strand text, new_start text, new_end text, new_strand text, new_locus_tag text)}

## Read in the weight file (typically many GBs) line by line ##
set number 0
set fp [open $file_to_sqlite_DB]
while 1 {
	
	gets $fp line
	## If eof is reached, then exit the while loop ##
	if {[eof $fp]} {break}
	
	## Otherwise process each line, extracting the entries ##
	set line [split [string trim $line] \t]

	## Ignore any blank lines
	if {[llength $line] == 0} {
		puts "ERROR: there are not 3 elements in this line:\t$line"
		continue
	}

	set old_id [lindex $line 12]

	## Don't input entries where there is no old ID - we will never query these ##
	if {[string length $old_id] == 0} {
		continue
	}

	set new_id [lindex $line 23]

	set accession [lindex $line 1]
	set old_start [lindex $line 3]
	set old_end [lindex $line 4]
	set old_strand [lindex $line 5]
	set new_start [lindex $line 14]
	set new_end [lindex $line 15]
	set new_strand [lindex $line 16]
	set new_locus_tag [lindex $line 20]

	#puts "$old_id\t$new_id\t$accession\t$old_start\t$old_end\t$old_strand\t$new_start\t$new_end\t$new_strand"

	## Insert the data into the sqlite database ##
	db1 eval {insert into t1 values($old_id,$new_id,$accession,$old_start,$old_end,$old_strand,$new_start,$new_end,$new_strand,$new_locus_tag)}
}

puts "Creating indexes for $file_to_sqlite_DB"
db1 eval {CREATE INDEX old_id_index ON t1 (old_id)}
db1 eval {CREATE INDEX new_id_index ON t1 (new_id)}

db1 close
