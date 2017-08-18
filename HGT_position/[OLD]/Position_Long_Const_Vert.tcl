#!/usr/local/bin/tclsh

source ~/Dropbox/Scripts/General_utils.tcl
package require sqlite3

proc geo_origin_positions {} {
	puts "Getting a table of the dnaA / dnaN gene positions to determine genome origin ..."

	set origin_genes_l [split [string trim [openfile /Users/aesin/Desktop/Geo_analysis/HGT_position/Geo_origin_genes.txt]] \n]

	global ori_rep_tbl {}
	global species_l {}

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

	puts "DONE: Table of origin positions"
	return
}

proc geo_assembled_only {} {

	if {[file exists /users/aesin/desktop/Geo_analysis/HGT_position/assembled_geo.tsv] == 1} {
		set assembled_species_l [split [string trim [openfile /users/aesin/desktop/Geo_analysis/HGT_position/assembled_geo.tsv]] \n]
	
	} else {
		cd /Users/aesin/Desktop/Geo_analysis/Geo_omes/Genome_processing/Geobac_genomes_assembled_genbank
		set non_conting_genomes [glob *gbk]
		regsub -all {_genomic.gbk} $non_conting_genomes {} non_contig_ncbi_ids

		set assembled_species_l {}
		foreach assembled_ncbi_id $non_contig_ncbi_ids {
			set assembled_species [db1 eval {select binomial from t1 where ncbi_id = $assembled_ncbi_id LIMIT 1}]
			lappend assembled_species_l $assembled_species
		}

		set out [open /users/aesin/desktop/Geo_analysis/HGT_position/assembled_geo.tsv w]
		puts $out [join $assembled_species_l \n]
		close $out
	}
	puts "DONE: Identified which species are fully assembled genomes"
	return $assembled_species_l
}

proc get_relative_position {start_pos species_start relevant_genome_length} {
	set ori_strand [lindex [split $species_start \t] 3]
	if {$ori_strand eq "+"} {
		set ori_start [lindex [split $species_start \t] 1]

		# Get the relative position #
		set relative_pos [expr $start_pos - $ori_start]
		# If it's negative then it lies upstream of the origin and we append it to the "end" #
		if {$relative_pos < 0} {
			set relative_pos [expr $relative_pos + $relevant_genome_length]
		}

	} else {
		set ori_start [lindex [split $species_start \t] 2]

		# Get the relative position #
		set relative_pos [expr $ori_start - $start_pos]
		# If it's negative then it lies upstream of the origin and we append it to the "end" #
		if {$relative_pos < 0} {
			set relative_pos [expr $relative_pos + $relevant_genome_length]
		}

	}
	
	## Get the fractional position ##
	set fraction_pos [tcl::mathfunc::roundto [expr double($relative_pos) / double($relevant_genome_length)] 5]
	return $fraction_pos
}

###########################################################################
## If linear representation, then linear == TRUE, if circular linear == FALSE

set penalty 5
set linear FALSE

## Open the ortholog database ##
cd /Users/aesin/Desktop/Geo_analysis/Geo_ortholog_nucl
sqlite3 db1 Geo_v_all_orth_location_NEW_db

###########################################################################
## Table of species with their dnaA and polB gene locus tags ##
geo_origin_positions

###########################################################################
## Which genomes are contig-based? ##
set assembled_species_l [geo_assembled_only]

###########################################################################
## Total genome lengths (adjusted for plasmids) ##
set genome_length_l [split [string trim [openfile  /Users/aesin/Desktop/Geo_analysis/Geo_omes/Genome_lengths.tsv]] \n]
puts "DONE: Table of genome lengths"

#######################################################################################
## Get transfer position data for long-distance and consistent sets at given penalty ##
#######################################################################################

set direct /users/aesin/Desktop/Mowgli

set receptor_genomes_dir $direct/Mowgli_outputs
set tip_key_dir $direct/OLD/Tip_keys

set receptor_table_file $receptor_genomes_dir/Per_penalty_tips/Per_penalty_tips_t$penalty\.tsv
set receptor_table_data [split [string trim [openfile $receptor_table_file]] \n]


set data_sets {long const vert all}

foreach set $data_sets {
	
	######################
	## Set output paths ##
	######################

	if {$linear == TRUE} {
		set output_dir /users/aesin/desktop/Geo_analysis/HGT_position/For_linear/Per_species_$set
	} else {
		set output_dir /users/aesin/desktop/Geo_analysis/HGT_position/For_circular/Per_species_$set
	}

	file mkdir $output_dir/Full_entries
	file mkdir $output_dir/Start
	file mkdir $output_dir/Relative_start

	############################
	## Set counters and lists ##
	############################

	set tried_to_map 0
	set failed_to_map 0
	set plasmid_gene 0

	set more_than_1_xfer_in_group 0
	set total_geo_ids 0

	set plasmid_entries {}
	set species_list {}

	set event_counter 1

	if {$set eq "const"} {

		cd $direct/Consistent_HGT_Vertical/Consistent_Full/Scenarios_1_2
		set HGT_intersect_file [lindex [glob true_HGT_intersect*$penalty\.tsv] 0]
		if {[string length $HGT_intersect_file] == 0} {puts "CAN'T FIND FILE"; exit}
		set events [lrange [split [string trim [openfile $HGT_intersect_file]] \n] 0 end]

	} elseif {$set eq "long"} {

		set long_distance_file $direct/Long_distance_HGT/Full/Scenarios_1_2/Donor_edges/T$penalty\_refined_receptor_donor_list.txt
		set events [lrange [split [string trim [openfile $long_distance_file]] \n] 0 end]
	
	} elseif {$set eq "vert"} {

		cd $direct/Consistent_HGT_Vertical/Consistent_Full/Scenarios_1_2/
		set vertical_file [lindex [glob true_vertical_intersect_t$penalty*t20.tsv] 0]
		set events [lrange [split [string trim [openfile $vertical_file]] \n] 0 end]

	} elseif {$set eq "all"} {

		set events [db1 eval {select DISTINCT group_number from t1 where binomial glob "Geobacillus*"}]
	}

	#########################################################################
	## For each assembled Geobacillus genome create a set of output tables ##
	#########################################################################
	foreach species $assembled_species_l {
		set $species\_tbl {}
		set $species\_start_tbl {}
		set $species\_rel_start_tbl {}
		lappend species_list $species
	}
	

	foreach event $events {

		## Process the input data "events" depending on the type of input ##
		if {$set eq "const" || $set eq "vert" || $set eq "all"} {
			
			set group $event

		} elseif {$set eq "long"} {

			set group [lindex [split [lindex [split $event \t] 1] " "] 0]
			set receptor_node [join [split [lindex [split $event \t] 0] " "] ","]
			set donor_node [join [lrange [split [lindex [split $event \t] 1] " "] 1 end] ","]
		
		}

		## If this is a transfer, we need to indentify a single unique event ##
		if {$set eq "long" || $set eq "const"} {

			## Open the tip key translator file ##
			set tip_transl_tbl [split [string trim [openfile $tip_key_dir/TIP_KEY_$group\.tsv]] \n]

			## Locate the corresponding entry in our table of receptor genomes ##
			set entry_in_recpt_tbl [lsearch -all -inline -glob $receptor_table_data $group\t*]
			if {[string length $entry_in_recpt_tbl] == 0} {
				puts "not found $group"
				break
			}

			if {$set eq "const"} {

				## If there is more than one HGT in this group, ignore it ##
				if {[llength $entry_in_recpt_tbl] > 1} {
					incr event_counter
					incr more_than_1_xfer_in_group
					continue
				}

			} elseif {$set eq "long"} {

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
					#incr num_xfers
				}

			}

			set entries_to_be_processed [split [lindex [split [lindex $entry_in_recpt_tbl 0] \t] end] " "]

		
		} elseif {$set eq "vert" || $set eq "all"} {

			set entries_to_be_processed [db1 eval {SELECT prot_id from t1 where group_number = $group AND binomial glob "Geobacillus*";}]
			set total_geo_ids [expr $total_geo_ids + [llength $entries_to_be_processed]]
		} 


		foreach entry $entries_to_be_processed {
			
			if {$set eq "long" || $set eq "const"} {
				set species_id [lsearch -all -inline -glob $tip_transl_tbl *\t$entry]

				if {[llength $species_id] != 1} {
					puts "Issue with finding the tip transation for group: $group and species tip: $entry"
					exit
				} else {
					set species_id [lindex [split [lindex $species_id 0] \t] 0] 
					set species_id [string range $species_id [string first "\}\{" $species_id]+2 end-1]
				}

			} elseif {$set eq "vert" || $set eq "all"} {
				set species_id $entry
			}
			
			incr tried_to_map
			set entry_l [split [db1 eval {SELECT * from t1 where prot_id = $species_id}] " "]
			if {[llength $entry_l] == 0} {
				incr failed_to_map
				continue
			}

			set species [lindex $entry_l 9]
			# Only continue if the genome is assembled #
			if {[lsearch $assembled_species_l $species] == -1} {
				continue
			} else {
				if {[lsearch $species_list $species] == -1} {
					lappend species_list $species
				}
			}

			################################################################################
			## Check the protein ID against the plasmid TAG list for the relevent species ##
			################################################################################

			## Open Plasmid tag file (if one exists for this species) ##
			set plasmid_tag_file /users/aesin/Desktop/Geo_analysis/Geo_omes/Plasmid_tags/$species\_plasmid_locus_tags.tsv
			if {[file exists $plasmid_tag_file] == 1} {
				set tag_list [split [string trim [openfile $plasmid_tag_file]] \n]
				set transfer_locus [lindex $entry_l 8]

				## If this gene DOES map to the plasmid, do NOT add it to the total output list ##
				if {[lsearch $tag_list $transfer_locus] != -1} {
					incr plasmid_gene
					lappend plasmid_entries $entry_l
					continue
				}
			}

			######################################################################################
			## Put the information (either full or just start position) into per-species tables ##
			######################################################################################

			## The strand will determine whether the true gene start position is the start or the end position as in the table. I.e. for complement strand, the "start" in the table is in fact the end, and vice versa ##
			## Start positions only ##
			set strand [lindex $entry_l 7]
			if {$strand eq "+"} {
				set start_pos [lindex $entry_l 5]
			} else {
				set start_pos [lindex $entry_l 6]
			}
	
			###########################################################################
			## Make the distance to the origin a relative value of the genome length ##
			###########################################################################

			## Relative start positions ##
			if {[regexp "Geobacillus" $species] == 1} {
				set relevant_genome_length [lindex [split [lsearch -inline -glob $genome_length_l $species\t*] \t] 1]

				## Get the origin position for this genome ##
				set species_start [lsearch -glob -inline $ori_rep_tbl $species\t*]

				## If the origin lies on the complement strand, we take the "end" as the start. We also change the upstream/downstream orientation. ##
				
				set fraction_pos [get_relative_position $start_pos $species_start $relevant_genome_length]
				
				## Add raw position to the raw start table ##
				lappend $species\_start_tbl $start_pos
				## Append the fractional position to the full entry table ##
				lappend $species\_tbl [join [concat $entry_l $fraction_pos] " "]
				## Add fractional position to the per-species table ##
				lappend $species\_rel_start_tbl $fraction_pos

				## If the position is in the first third or so of the genome, add 1 to have it rollover so that we can see the correct density across the origin ##
				if {$linear == TRUE && $fraction_pos < 0.3} {
					set rollover_fraction_pos [expr $fraction_pos + 1]
					lappend $species\_rel_start_tbl $rollover_fraction_pos	
				}
			}

		}

		#puts "Events sorted: $event_counter / [llength $events]"
		incr event_counter

	}

	set combined_rel_pos {}

	foreach species $species_list {
		set full_out_tbl [join [expr $$species\_tbl] \n]
		set start_out_tbl [join [expr $$species\_start_tbl] \n]
		
		if {$set eq "long" || $set eq "const"} {
			set output_name "$species\_HGT.tsv"
		} elseif {$set eq "vert"} {
			set output_name "$species\_vert.tsv"
		} elseif {$set eq "all"} {
			set output_name "$species\_all.tsv"
		}

		set out [open $output_dir/Full_entries/$output_name w]
		puts $out $full_out_tbl
		close $out

		set out [open $output_dir/Start/$output_name w]
		puts $out $start_out_tbl
		close $out

		if {[regexp "Geobacillus" $species] == 1} {
			set rel_start_out_tbl [join [expr $$species\_rel_start_tbl] \n]
			
			set out [open $output_dir/Relative_start/$output_name w]			
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
 
	puts "\n\n\#\#\#\#\#\# Set: $set \#\#\#\#\#\#"
	puts "Events processed: [expr $event_counter - 1]"
	if {$set eq "const"} {
		puts "Number of consistent transfer groups with > 1 transfer: $more_than_1_xfer_in_group"
	} elseif {$set eq "vert" || $set eq "all"} {
		puts "Total geo ids processed: $total_geo_ids"
	}
	puts "Genes we tried to map: $tried_to_map"
	puts "Genes that failed to map: $failed_to_map"
	puts "Genes mapping to plasmid: $plasmid_gene"

}

db1 close

