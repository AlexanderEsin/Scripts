#!/usr/local/bin/tclsh

source ~/Dropbox/Scripts/General_utils.tcl
package require sqlite3

###
set direct /users/aesin/desktop/Clean_Genomes
set db_name Feature_table_db
###

cd $direct

## Create the sqlite database within the Rbhh folder ##
sqlite3 db1 $db_name
db1 eval {PRAGMA main.page_size = 4096}
db1 eval {PRAGMA main.cache_size=10000}
db1 eval {PRAGMA main.locking_mode=EXCLUSIVE}
db1 eval {PRAGMA main.synchronous=NORMAL}
db1 eval {PRAGMA main.journal_mode=WAL}
db1 eval {create table t1(ncbi_id text, prot_id text, accession text, gene_start text, gene_end text, strand text, locus_tag text)}

## Go the the feature table folder ##

cd $direct/Prok_feature_tables
set feature_tables [glob *txt]
set num_tables [llength $feature_tables]
# set feature_tables [lrange $feature_tables 0 0]

set done_counter 0
set break_code 0

foreach table $feature_tables {

	set ncbi_id ""
	set ncbi_id [string range $table 0 [string first "_feature" $table]-1]
	puts $ncbi_id

	set feature_data [split [string trim [openfile $table]] \n]
	# set feature_data [lrange $feature_data 0 99]

	set line_index 0
	foreach line $feature_data {
		## Ignore comment lines ##
		if {[string index $line 0] ne "#"} {
			set entry [split $line \t]
			
			if {[lindex $entry 0] eq "gene" && [lindex $entry 1] eq "protein_coding"} {

				set accession ""
				set prot_id ""
				set gene_start ""
				set gene_end ""
				set strand ""
				set locus_tag ""
			
				set accession [lindex $entry 6]
				set gene_start [lindex $entry 7]
				set gene_end [lindex $entry 8]
				set strand [lindex $entry 9]

				set locus_tag [lindex $entry 16]

				############ Find the protein ID ############

				## Assume the CDS entry is the next line ##
				set CDS_entry [split [lindex $feature_data $line_index+1] \t]
				set CDS_locus_tag [lindex $CDS_entry 16]

				## Let's check that the next line is 1. a CDS entry and 2. has the same locus tag ##
				if {$locus_tag eq $CDS_locus_tag && [lindex $CDS_entry 0] eq "CDS"} {
					set prot_id [lindex $CDS_entry 10]

				## If not, let's try and find a matching locus tag CDS entry ##
				} else {
					set CDS_entry [split [lsearch -inline -glob $feature_data CDS*$locus_tag*] \t]
					set prot_id [lindex $CDS_entry 10]

					## If we can't find a matching CDS locus tag, go back to assuming it's simply the next line on, just make sure it's still a CDS line ##
					if {[string length $prot_id] == 0} {
						set CDS_entry [split [lindex $feature_data $line_index+1] \t]

						if {[lindex $CDS_entry 0] eq "CDS"} {
							set prot_id [lindex $CDS_entry 10]
							## If the gene locus_tag cannot be matched to the CDS locus tag, we will take the protein locus tag instead ##
							set locus_tag [lindex $CDS_entry 16]

						} else {
							puts stderr "Bollocks"
							puts stderr "$ncbi_id\t$prot_id\t$accession\t$gene_start\t$gene_end\t$strand\t$locus_tag"
							exit
						}
						
					}
				}

				#puts "$prot_id\t$gene_start\t$gene_end"
				
				##  If any of the necessary information is lacking - exit ##
				if {$ncbi_id eq "" || $prot_id eq "" || $accession eq "" || $gene_start eq "" || $gene_end eq "" || $strand eq ""|| $locus_tag eq ""} {
					puts stderr "ERROR - essential data missing"
					puts stderr "Locus tag: $locus_tag"
					puts stderr "$ncbi_id\t$prot_id\t$accession\t$gene_start\t$gene_end\t$strand\t$locus_tag"
					set break_code 1
					break
				} else {
					#puts "$ncbi_id\t$prot_id\t$accession\t$gene_start\t$gene_end\t$strand"
					db1 eval {insert into t1 values($ncbi_id,$prot_id,$accession,$gene_start,$gene_end,$strand,$locus_tag)}
				}
			}
		}
		incr line_index
	}

	if {$break_code == 1} {
		break
	}
	puts "Entered into DB: $done_counter / $num_tables"
	incr done_counter
}

if {$break_code != 1} {
	puts "Creating indexes for $db_name"
	db1 eval {CREATE INDEX old_id_index ON t1 (ncbi_id)}
	db1 eval {CREATE INDEX new_id_index ON t1 (prot_id)}
	db1 eval {CREATE INDEX locus_tag_index ON t1 (locus_tag)}

	db1 close
}


