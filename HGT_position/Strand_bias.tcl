#!/usr/local/bin/tclsh

source ~/Dropbox/Scripts/General_utils.tcl
source ~/Dropbox/Scripts/Procs/gc_content.tcl

package require sqlite3

###########################################################################
set penalty 6

## Table of species with their dnaA and polB gene locus tags ##
set origin_genes_l [split [string trim [openfile /Users/aesin/Desktop/Geo_analysis/HGT_position/Geo_origin_genes.txt]] \n]

## Table of groups with HGT and the receptor edges ##
set receptor_group_donor_tbl [split [string trim [openfile /Users/aesin/Desktop/Mowgli/Long_distance_HGT/Full/Scenarios_1_2/Donor_edges/T6_refined_receptor_donor_list.txt]] \n]

## Open the ortholog database ##
cd /Users/aesin/Desktop/Geo_analysis/Geo_ortholog_nucl
sqlite3 db1 Geo_v_all_orth_location_NEW_db


###########################################################################
## Identify on which strand the dnaA/polB genes lie - this is the origin of replication, and our strand data is relative to it. I.e. if the gene is on the same strand as the origin cassette - we give it a 1, if on a different strand, we give it a 0. ##

set reference_strand_tbl {}
set species_l {}

foreach entry $origin_genes_l {
	set entry_l [split $entry \t]
	set species [lindex $entry_l 0]
	set dnaA_gene [lindex $entry_l 1]

	set dnaA_strand [db1 eval {select strand from t1 where locus_tag = $dnaA_gene}]
	if {[string length $dnaA_strand] == 0} {
		set polB_gene [lindex $entry_l 2]
		set polB_strand [db1 eval {select strand from t1 where locus_tag = $polB_gene}]
		set initiation_strand $polB_strand
	} else {
		set initiation_strand $dnaA_strand
	}

	lappend reference_strand_tbl "$species\t$initiation_strand"
	lappend species_l $species
}

###########################################################################
### Mowgli can't have numbers at the end of the species name (e.g. ...thermodenitrificans_NG80_2), so they were converted (e.g. ...thermodenitrificans_NG80_2A). This part converts the mowgli-parsable names back into the original species names that are present in the database ###

## Table translating nodes (to which HGT is predicted) and all the subtended Anoxy/Geobacillus tips under that node ##
set node_tip_transl_l [lrange [split [string trim [openfile /Users/aesin/Desktop/Geo_analysis/HGT_position/Geo_nodes_tips.txt]] \n] 1 end]
set node_tip_new_l $node_tip_transl_l

foreach species $species_l {
	set pat \{\(?:$species\)+\[A-Z\]*\}
	set species_in_tree [regexp -inline [expr $pat] $node_tip_transl_l]
	puts "$species_in_tree\t$species"
	regsub -all $species_in_tree $node_tip_new_l $species node_tip_new_l
}

###########################################################################
## Get a list of all the groups in the database ##
# set groups [lsort -dictionary [db1 eval {SELECT DISTINCT group_number from t1}]]

set strand_out_tbl {}
set break_point 0

## Counters ##
set total_genes 0
set mapped 0
set not_mapped 0
set same_strand 0
set diff_strand 0

#set receptor_group_donor_tbl [lrange $receptor_group_donor_tbl 0 5]
foreach hgt_event $receptor_group_donor_tbl {
	set receptor_group_donor [split $hgt_event \t]
	set group [lindex [split [lindex $receptor_group_donor 1] " "] 0]
	set receptor_node [lindex [split [lindex $receptor_group_donor 0] " "] 1]

	set line_output $group

	set subtended_species [split [lindex [split [lsearch -inline -glob $node_tip_new_l $receptor_node\t*] \t] 1] " "]

	foreach species $species_l {

		if {[lsearch $subtended_species $species] == -1} {lappend line_output {.}; continue}

		incr total_genes

		set entry [db1 eval {SELECT * from t1 where group_number = $group AND binomial = $species}]

		if {[string length $entry] == 0} {lappend line_output {.}; incr not_mapped; continue}

		incr mapped

		set entry_l [split $entry " "]
		set strand [lindex $entry_l 7]

		set reference_strand [lindex [split [lsearch -glob -inline $reference_strand_tbl $species\t*] \t] 1]
		if {$strand == $reference_strand} {
			lappend line_output 1
			incr same_strand
		} else {
			lappend line_output 0
			incr diff_strand
		}

	}
	if {$break_point == 1} {break}
	set line_output_str [join $line_output \t]
	lappend strand_out_tbl $line_output_str
}

puts "Total HGT genes: $total_genes\nTotal mapped: $mapped\nNot mapped: $not_mapped\n\nSame strand as dnaA/polB: $same_strand\nDifferent strand: $diff_strand"

#puts [join $strand_out_tbl \n]
#db1 close