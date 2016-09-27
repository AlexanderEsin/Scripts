#!/usr/local/bin/tclsh
source ~/Dropbox/Scripts/General_utils.tcl
package require sqlite3

## Make an sqlite db to store the cluster --> annotation info ##
## The .annotations.tsv file has the following (unwritten) header: TaxonomicLevel | GroupName | ProteinCount | SpeciesCount | COGFunctionalCategory | ConsensusFunctionalDescription ##

set direct /users/aesin/desktop
set working_directory $direct/Geo_analysis/Geo_ortholog_nucl/EggNOG

cd $working_directory

## Prepare the sqlite3 database ##
sqlite3 db1 eggNOG_annot_db
db1 eval {PRAGMA main.page_size = 4096}
db1 eval {PRAGMA main.cache_size=10000}
db1 eval {PRAGMA main.locking_mode=EXCLUSIVE}
db1 eval {PRAGMA main.synchronous=NORMAL}
db1 eval {PRAGMA main.journal_mode=WAL}
db1 eval {create table t1(group_id text, protein_count text, species_count text, COGcategory text, function blob)}


set annot_file bactNOG.annotations.tsv
set annot_entries [split [string trim [openfile $annot_file]] \n]

progress_init [llength $annot_entries]
## There are 144498 entries ##

set counter 0

foreach entry $annot_entries {
	set columns [split $entry \t]
	set group [lindex $columns 1]
	set protein_count [lindex $columns 2]
	set species_count [lindex $columns 3]
	set COGcategory [lindex $columns 4]
	set function [lindex $columns 5]

	db1 eval {insert into t1 values($group,$protein_count,$species_count,$COGcategory,$function)}
	incr counter
	progress_tick $counter
}

db1 eval {CREATE INDEX old_id_index ON t1 (group_id)}

db1 close