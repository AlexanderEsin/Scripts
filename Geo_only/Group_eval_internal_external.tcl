#!/usr/local/bin/tclsh

#####################################################################################################################################
# This script:																														#
# 		1) makes an sqlite database of all the Rbbh_weight files that are used as input for clustering								#
#		2) compares the internal average evalues for a particular family (e.g. geo_only groups) and 								#
#		   identifies any external evalues to other families (i.e. show any potential clustering artefacts) 						#
#																																	#
#####################################################################################################################################

source ~/Dropbox/Scripts/General_utils.tcl
package require sqlite3
###########################################################################
## Switch to make the sqlite databases for the rbbh master weight files ##
set make_sqlite_db_switch 0

## Set the evalue at which to query the evalue hits (lowest stringency is most expansive and probably makes most sense) ##
set query_eval -10

## Set ival at which the analysis was done ##
set ival 2.0

## Path to master directory containing analysis ##
set direct /users/aesin/desktop/Geo_analysis/Geo_v_all

## Path to directory for the Geo_only group analysis ##
set geo_only_dir $direct/$ival/Geo_only_groups

###########################################################################
## Make the sqlite DBs for each of the Rbbh_Master_weight files ##

if {$make_sqlite_db_switch == 1} {
	cd $direct

	set rbbh_dirs [glob -type d Rbbh*]

	foreach dir $rbbh_dirs {
		cd $direct/$dir

		puts "Adding to DB for $dir ...."

		## Trim the folder name to extract the numeric eval ##
		set eval [string range $dir 9 end]

		## Create the sqlite database within the Rbhh folder ##
		set db_name Master_weight_$eval\_db
		sqlite3 db1 $db_name
		db1 eval {PRAGMA main.page_size = 4096}
		db1 eval {PRAGMA main.cache_size=10000}
		db1 eval {PRAGMA main.locking_mode=EXCLUSIVE}
		db1 eval {PRAGMA main.synchronous=NORMAL}
		db1 eval {PRAGMA main.journal_mode=WAL}
		db1 eval {create table t1(Gene_A blob, Gene_B blob, Eval blob)}

		set weight_file [lindex [glob *weight.txt] 0]

		if {[llength $weight_file] != 1} {
			puts "ERROR: Cannot find weight file!"
			exit 2
		}

		## Read in the weight file (typically many GBs) line by line ##
		set number 0
		set fp [open $weight_file]
		while 1 {
			
			gets $fp line
			## If eof is reached, then exit the while loop ##
			if {[eof $fp]} {break}
			
			## Otherwise process each line, extracting the entries ##
			set line [split [string trim $line] \t]
			if {[llength $line] != 3} {
				## Sometimes there seems to be a blank line, perhaps an issue with data insertion at the rbbh_identification level. These are ignored ##
				if {$line != ""} {
					puts "ERROR: there are not 3 elements in this line:\t$line"
				}
				continue
			}

			set gene_1 [lindex $line 0]
			set gene_2 [lindex $line 1]
			set eval_temp [lindex $line 2]

			## Insert the data into the sqlite database ##
			db1 eval {insert into t1 values($gene_1,$gene_2,$eval_temp)}
		}

		## Create an index on both the gene columns, even though we only query the first column (Gene_A) downstream ##
		puts "Creating indexes for $dir"
		db1 eval {CREATE INDEX gene_A_index ON t1 (Gene_A)}
		db1 eval {CREATE INDEX gene_B_index ON t1 (Gene_B)}

		db1 close
	}
}
###########################################################################
## Perform the evalue analysis ##

cd $geo_only_dir/Group_fastas
set geo_only_fastas [glob *faa]

set rbbh_query_database $direct/Rbbh_all_$query_eval/Master_weight_$query_eval\_db
sqlite3 db1 $rbbh_query_database

set groups_query_database $direct/$ival/$query_eval/Groups_$query_eval\_db
sqlite3 db2 $groups_query_database

set output_table {}
lappend output_table "Group_no\tNumber of Geobac genes\tMax possible number of internal hits\tNumber of internal hits\tAverage internal evalue\tNumber of external hits\tAverage external evalue\tRelated groups\tExternal IDs"
set i 1

#set geo_only_fastas [lrange $geo_only_fastas 0 0]

foreach fasta $geo_only_fastas {
	set fasta_group_no [string range $fasta 0 end-4]
	set data [string trim [openfile $fasta]]
	set genes [split_genes_fast $data]

	## Make a list of all the gene ids in this group ##
	set geo_only_ids {}
	foreach gene $genes {
		set gene_id [string trim [string range $gene 1 [string first " " $gene]]]
		lappend geo_only_ids $gene_id
	}

	set geo_gene_no [llength $genes]
	set max_possible_internal_hits [expr [llength $geo_only_ids] * ([llength $geo_only_ids] - 1)]

	## Query the database for each gene id to find all hits - then filter those against ids in this family separately ##
	set internal_hits_evals {}
	set external_hits_evals {}
	set external_hits_ids {}
	set total_hits 0

	foreach id $geo_only_ids {

		## Find all the ids and respective evalues that were hit by this id ##
		set hits [db1 eval {SELECT (Gene_B||" "|| Eval) AS hit FROM t1 WHERE Gene_A = $id}]
		set total_hits [expr $total_hits + [llength $hits]]

		## The resultant list from above looks like this {Gene_B1 Eval1} {Gene_B2 Eval2} ... {Gene_Bn Evaln} ##

		foreach hit $hits {
			set gene_b [lindex [split $hit " "] 0]
			set hit_eval [lindex [split $hit " "] 1]
			## Here we split each hit, and then compare the Gene_B entry to the list of genes in this gene family. If the hit is to a member of the gene family, then it's counted as internal. Otherwise counted as external ##
			set internal_counter 0

			if {[lsearch $geo_only_ids $gene_b] > -1} {
				lappend internal_hits_evals $hit_eval
				incr internal_counter
			} else {
				lappend external_hits_ids $gene_b
				lappend external_hits_evals $hit_eval
			}
		}
		## The number of internal hits should not exceed the number of ids in the family - 1 ##
		if {$internal_counter > [expr [llength $geo_only_ids] -1]} {
			puts "ERROR: more internal hits than should be possible\nFasta: $fasta\tID:$id\tInternal hits: $internal_counter\tMax possible: [expr [llength $geo_only_ids] -1]"
		}
	}

	set no_internal_hits [llength $internal_hits_evals]
	set no_external_hits [llength $external_hits_evals]

	if {$total_hits == 0} {
		puts "WHAT?"
		exit 2
	}

	## Calculate the average internal evalue ##

	set log_eval_l {}
	foreach int_eval $internal_hits_evals {
		if {$int_eval == "0.0"} {
			set int_eval 1e-200
		}
		set log_eval [expr log($int_eval)]
		lappend log_eval_l $log_eval
	}
	set total_int_eval [ladd $log_eval_l]
	set av_int_eval [format %3.3g [expr exp($total_int_eval/$no_internal_hits)]]

	## If there are any external hits calculate their av eval and also join together a list of external IDs ##
	set av_ext_eval "-"
	set external_id_str "-"
	set related_groups_str "-"
	if {$no_external_hits > 0} {
		set log_eval_l {}
		foreach ext_eval $external_hits_evals {
			if {$ext_eval == "0.0"} {
				set ext_eval 1e-200
			}
			set log_eval [expr log($ext_eval)]
			lappend log_eval_l $log_eval
		}
		set total_ext_eval [ladd $log_eval_l]
		set av_ext_eval [format %3.3g [expr exp($total_ext_eval/$no_external_hits)]]

		set external_hits_ids [lsort -unique $external_hits_ids]
		set external_id_str [join $external_hits_ids " "]

		## For each external hit, see if they belong to any other group(s) ##
		set related_groups_l {}
		foreach external_id $external_hits_ids {
			set external_id_group ""
			set glob_pat "$external_id*"

			db2 eval {SELECT group_no FROM t1 WHERE genes GLOB $glob_pat} {
			    set external_id_group "$group_no"
			}

			if {$external_id_group != ""} {
				lappend related_groups_l $external_id_group
			}
		}
		if {[llength $related_groups_l] != 0} {
			set related_groups_l [lsort -unique $related_groups_l]
			set related_groups_str [string trim [join $related_groups_l " "]]
		}
	}
	
	## Add a line into the output table with the group stats ##
	lappend output_table "$fasta_group_no\t$geo_gene_no\t$max_possible_internal_hits\t$no_internal_hits\t$av_int_eval\t$no_external_hits\t$av_ext_eval\t$related_groups_str\t$external_id_str"

	puts "Groups done ... $i / [llength $geo_only_fastas]"
	incr i
}

set output_table_str [join $output_table \n]
set out [open $geo_only_dir/Geo_only_eval_stats.tsv w]
puts $out $output_table_str
close $out

db1 close
db2 close
###########################################################################



