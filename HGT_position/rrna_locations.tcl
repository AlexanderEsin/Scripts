#!/usr/local/bin/tclsh

source ~/Dropbox/Scripts/General_utils.tcl
source ~/Dropbox/Scripts/Procs/gc_content.tcl
package require sqlite3

###########################################################################
set direct /users/aesin/Desktop/Geo_analysis/Geo_omes/Assembled_geo_no_plasmids/Assembled_no_plasmid_gbks
set output_direct /users/aesin/Desktop/Geo_analysis/HGT_position/Ribosmal_rna_locations
file mkdir $output_direct/Full_tables
file mkdir $output_direct/Relative_start

cd /Users/aesin/Desktop/Geo_analysis/Geo_ortholog_nucl
sqlite3 db1 Geo_v_all_orth_location_NEW_db

###############################################################
## Table of species with their dnaA and polB gene locus tags ##
###############################################################

puts "Getting a table of the dnaA / dnaN gene positions to determine genome origin ..."

set origin_genes_l [split [string trim [openfile /Users/aesin/Desktop/Geo_analysis/HGT_position/Geo_origin_genes.txt]] \n]

set ori_rep_tbl {}
set species_l {}

foreach entry $origin_genes_l {
	set entry_l [split $entry \t]
	set species [lindex $entry_l 0]

	set dnaA_gene [lindex $entry_l 1]
	set dnaA_start_end [db1 eval {select gene_start, gene_end, strand from t1 where locus_tag = $dnaA_gene}]

	## If there's no entry for dnaA, resort to polB (dnaN) ##
	if {[string length $dnaA_start_end] == 0} {
		set polB_gene [lindex $entry_l 2]
		set polB_start_end [db1 eval {select gene_start, gene_end, strand from t1 where locus_tag = $polB_gene}]

		set polb_position_zero	[lindex $polB_start_end 0]
		set polb_position_last	[lindex $polB_start_end 1]
		set polb_strand			[lindex $polB_start_end 2]
		lappend ori_rep_tbl "$species\t$polb_position_zero\t$polb_position_last\t$polb_strand"
	} else {
		set dnaa_position_zero	[lindex $dnaA_start_end 0]
		set dnaa_position_last	[lindex $dnaA_start_end 1]
		set dnaa_strand			[lindex $dnaA_start_end 2]
		lappend ori_rep_tbl "$species\t$dnaa_position_zero\t$dnaa_position_last\t$dnaa_strand"
	}

	lappend species_l $species
}

set out [open /Users/aesin/Desktop/Geo_analysis/HGT_position/Geo_origin_positions.tsv w]
puts $out [join $ori_rep_tbl \n]
close $out
puts [join $ori_rep_tbl \n]
puts "DONE: Table of origin positions\n"


###########################################################################
set genome_length_l [split [string trim [openfile  /Users/aesin/Desktop/Geo_analysis/Geo_omes/Genome_lengths.tsv]] \n]
puts "DONE: Table of genome lengths\n"

###########################################################################
cd $direct
set no_plasmid_gbks [glob *gbk]

foreach gbk_file $no_plasmid_gbks {
	
	set ncbi [string range $gbk_file 0 [string last "_" $gbk_file]-1]
	set binomial [db1 eval {select binomial from t1 where ncbi_id = $ncbi LIMIT 1}]

	set relevant_genome_length [lindex [split [lsearch -inline -glob $genome_length_l $binomial\t*] \t] 1]

	set gbk_data [string trim [openfile $gbk_file]]
	set rrna_hits [regexp -all -inline {(?:rRNA)+?\s{3}.+?(?:product)+?.+?\n+?} $gbk_data]

	set rrna_output_table {}
	set relative_start_locations {}

	set species_start [lsearch -glob -inline $ori_rep_tbl $binomial\t*]
	set ori_strand [lindex [split $species_start \t] 3]
	if {$ori_strand eq "+"} {
		set ori_start [lindex [split $species_start \t] 1]
	} else {
		set ori_start [lindex [split $species_start \t] 2]
	}

	foreach rrna $rrna_hits {
		set rrna [string trim $rrna]

		## Get strand, start, and end of the gene ##
		set location [string trim [string range $rrna [string first " " $rrna] [string first \n $rrna]]]
		if {[regexp "complement" $location] != 1} {
			set strand "+"
			set coordinates [split [string trim $location] {..}]
			set coordinates [lremove_empty $coordinates]

			# Get the relative position #
			set start_pos [lindex $coordinates 0]
			set relative_pos [expr $start_pos - $ori_start]
			if {$relative_pos < 0} {
				set relative_pos [expr $relative_pos + $relevant_genome_length]
			}

			set fraction_pos [tcl::mathfunc::roundto [expr double($relative_pos) / double($relevant_genome_length)] 5]

		} else {
			set strand "-"
			set coordinates [split [string range $location [string first \( $location]+1 [string last \) $location]-1] {..}]
			set coordinates [lremove_empty $coordinates]

			set start_pos [lindex $coordinates 1]
			set relative_pos [expr $ori_start - $start_pos]
			# If it's negative then it lies upstream of the origin and we append it to the "end" #
			if {$relative_pos < 0} {
				set relative_pos [expr $relative_pos + $relevant_genome_length]
			}
			set fraction_pos [tcl::mathfunc::roundto [expr double($relative_pos) / double($relevant_genome_length)] 5]
		}

		set gene_start [lindex $coordinates 0]
		set gene_end [lindex $coordinates 1]

		lappend relative_start_locations $fraction_pos

		#puts "$strand\t$gene_start\t$gene_end"

		## Get name of the product ##
		if {[regexp -line "(?:product)+.+" $rrna product] == 0} {
			puts "No product identified\n$rrna"
		} else {
			set product [string range $product [string first \" $product]+1 [string last \" $product]-1]
		}
		
		## Get the new locus ID tag ##
		if {[regexp -line "(?:/locus_tag)+.+" $rrna locus_tag] == 0} {
			puts "No locus_tag identified\n$rrna"
		} else {
			set locus_tag [string range $locus_tag [string first \" $locus_tag]+1 [string last \" $locus_tag]-1]
		}

		lappend rrna_output_table "$binomial\t$gene_start\t$gene_end\t$strand\t$fraction_pos\t$locus_tag\t$product"
	}

	set out [open $output_direct/Full_tables/$binomial\_rrna_locations.tsv w]
	puts $out [join $rrna_output_table \n]
	close $out

	set out [open $output_direct/Relative_start/$binomial\_rrna_rel_start.tsv w]
	puts $out [join $relative_start_locations \n]
	close $out

}

db1 close