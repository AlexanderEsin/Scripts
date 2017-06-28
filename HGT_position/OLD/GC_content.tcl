#!/usr/local/bin/tclsh

source ~/Dropbox/Scripts/General_utils.tcl
package require sqlite3

###########################################################################
set direct /users/aesin/Desktop/Mowgli

set receptor_genomes_dir $direct/Mowgli_outputs
set tip_key_dir $direct/OLD/Tip_keys

set output_dir /users/aesin/desktop/Geo_analysis/HGT_position/Per_species_long
file mkdir $output_dir/Full_entries
file mkdir $output_dir/Start
file mkdir $output_dir/Relative_start

## Open the ortholog database ##
cd /Users/aesin/Desktop/Geo_analysis/Geo_ortholog_nucl
sqlite3 db1 Geo_v_all_orth_location_NEW_db

###########################################################################

set long_dist_xfers [lrange [split [string trim [openfile $direct/Long_distance_HGT/Full/Scenarios_1_2/Donor_edges/T5_refined_receptor_donor_list.txt]] \n] 0 end]

set receptor_table_file $receptor_genomes_dir/Test_output_penalty5.tsv
set receptor_table_data [split [string trim [openfile $receptor_table_file]] \n]

###########################################################################

## Table of species with their dnaA and polB gene locus tags ##
set origin_genes_l [split [string trim [openfile /Users/aesin/Desktop/Geo_analysis/HGT_position/Geo_origin_genes.txt]] \n]

set ori_rep_tbl {}
set species_l {}

foreach entry $origin_genes_l {
	set entry_l [split $entry \t]
	set species [lindex $entry_l 0]
	set dnaA_gene [lindex $entry_l 1]

	set dnaA_start [db1 eval {select gene_start from t1 where locus_tag = $dnaA_gene}]
	if {[string length $dnaA_start] == 0} {
		set polB_gene [lindex $entry_l 2]
		set polB_start [db1 eval {select gene_start from t1 where locus_tag = $polB_gene}]
		set position_zero $polB_start
	} else {
		set position_zero $dnaA_start
	}

	lappend ori_rep_tbl "$species\t$position_zero"
	lappend species_l $species
}

###########################################################################
## Total genome lengths (adjusted for plasmids) ##
set genome_length_l [split [string trim [openfile  /Users/aesin/Desktop/Geo_analysis/Geo_omes/Genome_lengths.tsv]] \n]

###########################################################################
## Which genomes are contig-based? ##
cd /Users/aesin/Desktop/Geo_analysis/Geo_omes/Genome_processing/Geobac_genomes_assembled_genbank
set non_conting_genomes [glob *gbk]
regsub -all {_genomic.gbk} $non_conting_genomes {} non_contig_ncbi_ids

set assembled_species_l {}
foreach assembled_ncbi_id $non_contig_ncbi_ids {
	set assembled_species [db1 eval {select binomial from t1 where ncbi_id = $assembled_ncbi_id LIMIT 1}]
	lappend assembled_species_l $assembled_species
}
###########################################################################

set num_xfers 0
set tried_to_map 0
set failed_to_map 0
set plasmid_xfer_gene 0

set plasmid_entries {}
set species_list {}

set xfer_counter 1

foreach transfer $long_dist_xfers {
	set group [lindex [split [lindex [split $transfer \t] 1] " "] 0]
	set receptor_node [join [split [lindex [split $transfer \t] 0] " "] ","]
	set donor_node [join [lrange [split [lindex [split $transfer \t] 1] " "] 1 end] ","]
	

	## Locate the corresponding entry in our table of receptor genomes ##
	set entry_in_recpt_tbl [lsearch -all -inline -glob $receptor_table_data $group\t*]
	if {[string length $entry_in_recpt_tbl] == 0} {
		puts "not found $group"
		break
	}

	## Ensure we find the specific transfer event in that group (in case there was more than one transfer into Geobacillus in that group), using the receptor & donor nodes ##
	if {[llength $entry_in_recpt_tbl] > 1} {
		#puts $entry_in_recpt_tbl
		set entry_in_recpt_tbl [lsearch -all -inline -glob $entry_in_recpt_tbl *$donor_node*]
		if {[llength $entry_in_recpt_tbl] > 1} {
			set entry_in_recpt_tbl [lsearch -all -inline -glob $entry_in_recpt_tbl *$receptor_node*]
		}
	}

	## If there is still ambiguity, i.e. > 1 event within the same group, same donor clade and same receptor node, break ##
	if {[llength $entry_in_recpt_tbl] > 1} {
		puts "There is still more than one entry that matches for $group group"
		break
	} else {
		set entry_in_recpt_tbl [lindex $entry_in_recpt_tbl 0]
		incr num_xfers
	}

	#########################################################################
	## For each transfer event establish the recipient Geobacillus species ##
	#########################################################################

	set tip_transl_tbl [split [string trim [openfile $tip_key_dir/TIP_KEY_$group\.tsv]] \n]

	set species_w_xfer [split [lindex [split $entry_in_recpt_tbl \t] end] " "]
	foreach species $species_w_xfer {
		set species_id [lsearch -all -inline -glob $tip_transl_tbl *\t$species]

		if {[llength $species_id] != 1} {
			puts "Issue with finding the tip transation for group: $group and species tip: $species"
			break
		} else {
			set species_id [lindex [split [lindex $species_id 0] \t] 0] 
			set species_id [string range $species_id [string first "\}\{" $species_id]+2 end-1]
		}

		incr tried_to_map
		set entry_l [split [db1 eval {SELECT * from t1 where prot_id = $species_id}] " "]
		if {[llength $entry_l] == 0} {
			incr failed_to_map
			continue
		}

		set species [lindex $entry_l 9]

		################################################################################
		## Check the protein ID against the plasmid TAG list for the relevent species ##
		################################################################################

		## Open Plasmid tag file (if one exists for this species) ##
		set plasmid_tag_file /users/aesin/Desktop/Geo_analysis/Geo_omes/Plasmid_tags/$species\_plasmid_locus_tags.tsv
		if {[file exists $plasmid_tag_file] == 1} {
			set tag_list [split [string trim [openfile $plasmid_tag_file]] \n]
			set transfer_locus [lindex $entry_l 8]

			## If this gene DOES map to the plasmid, do NOT add it to the total output list ##
			if {[lsearch $tag_list $transfer_locus] == 1} {
				incr plasmid_xfer_gene
				lappend plasmid_entries $entry_l
				continue
			}
		}

		######################################################################################
		## Put the information (either full or just start position) into per-species tables ##
		######################################################################################

		## Full table ##
		if {[info exists $species\_tbl] == 0} {
			set $species\_tbl {}
			lappend species_list $species
		} 
		lappend $species\_tbl [join $entry_l " "]

		## Start positions only ##
		set start_pos [lindex $entry_l 5]
		###### TEMPORARY ######
		## Due to issues with the ortholog table, remove any position entry that is not an integer ##
		if {[string is integer $start_pos] == 0} {
			puts "not integer: $start_pos"
			continue
		}
		if {[info exists $species\_start_tbl] == 0} {
			set $species\_start_tbl {}
		} 
		lappend $species\_start_tbl $start_pos

		###########################################################################
		## Make the distance to the origin a relative value of the genome length ##
		###########################################################################

		## Relative start positions ##
		if {[regexp "Geobacillus" $species] == 1} {
			set species_start [lsearch -glob -inline $ori_rep_tbl $species\t*]
			set ori_start [lindex [split $species_start \t] 1]
			set relative_pos [expr abs($start_pos - $ori_start)]

			## Get the fractional position ##
			set relevant_genome_length [lindex [split [lsearch -inline -glob $genome_length_l $species\t*] \t] 1]
			set fraction_pos [tcl::mathfunc::roundto [expr double($relative_pos) / double($relevant_genome_length)] 5]
			# if {$fraction_pos >= 0.5} {
			# 	set fraction_pos [tcl::mathfunc::roundto [expr abs($fraction_pos - double(1))] 5]
			# }
			if {$fraction_pos < 0.3} {

			}

			if {[info exists $species\_rel_start_tbl] == 0} {
				set $species\_rel_start_tbl {}
			} 
			lappend $species\_rel_start_tbl $fraction_pos

			## If the position is in the first third or so of the genome, add 1 to have it rollover so that we can see the correct density across the origin ##
			if {$fraction_pos < 0.3} {
				set rollover_fraction_pos [expr $fraction_pos + 1]
				lappend $species\_rel_start_tbl $rollover_fraction_pos	
			}
		}

	}

	puts "Transfer events sorted: $xfer_counter / [llength $long_dist_xfers]"
	incr xfer_counter
}

set combined_rel_pos {}
foreach species $species_list {
	set full_out_tbl [join [expr $$species\_tbl] \n]
	set start_out_tbl [join [expr $$species\_start_tbl] \n]
	

	set out [open $output_dir/Full_entries/$species\_HGT.tsv w]
	puts $out $full_out_tbl
	close $out

	set out [open $output_dir/Start/$species\_HGT.tsv w]
	puts $out $start_out_tbl
	close $out

	if {[regexp "Geobacillus" $species] == 1} {
		set rel_start_out_tbl [join [expr $$species\_rel_start_tbl] \n]
		set out [open $output_dir/Relative_start/$species\_HGT.tsv w]
		puts $out $rel_start_out_tbl
		close $out

		if {[lsearch $assembled_species_l $species] != -1} {
			lappend combined_rel_pos $rel_start_out_tbl
		}
	}	
}

set out [open $output_dir/combined_rel_pos_assembled.tsv w]
puts $out [join $combined_rel_pos \n]
close $out

db1 close

puts "\nXfers processed: $num_xfers"
puts "Genes we tried to map: $tried_to_map"
puts "Genes that failed to map: $failed_to_map"
puts "Genes mapping to plasmid: $plasmid_xfer_gene"
puts "Plasmid entries:\n[join $plasmid_entries \n]"