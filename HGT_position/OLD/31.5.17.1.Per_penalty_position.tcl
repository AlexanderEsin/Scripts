#!/usr/local/bin/tclsh

source ~/Documents/Scripts/General_utils.tcl
package require sqlite3

# Retrieve gene start/end positions and strand from sqlite DB #
proc GetGenePosByLocusTag {locus_tag {sqdb db1} {sqtbl t1}} {
	set gene_pos_data [$sqdb eval "SELECT gene_start, gene_end, strand FROM $sqtbl WHERE locus_tag = :locus_tag"]
	return $gene_pos_data
}

# For each species, identify the "Origin of Replication" - this is assumed equivalent to the positions of the DnaA gene (or PolB if DnaA absent) #
proc GetGeoOriPositions {} {
	
	# Input table format (by row): species_binomial [tab] DnaA_gene_locus_tag [tab] PolB_gene_locus_tag
	set origin_genes_l [split [string trim [openfile /Users/aesin/Desktop/Geo_analysis/HGT_position/Geo_origin_genes.txt]] \n]

	# Output table format: species_binomial [tab] gene_start_pos [tab] gene_end_pos [tab] strand
	set output_tbl {}

	foreach row $origin_genes_l {
		set columns	[split $row \t]
		set species	[lindex $columns 0]
		set dnaA_locus_tag	[lindex $columns 1]

		# Get position data of dnaA gene by locus_tag
		set dnaA_pos		[GetGenePosByLocusTag $dnaA_locus_tag]
		# If there's no position data for dnaA, try polB (dnaN)
		if {[string length $dnaA_pos] == 0} {
			set polB_locus_tag	[lindex $columns 2]
			set polB_pos		[GetGenePosByLocusTag $polB_locus_tag]

			if {[string length $polB_pos] == 0} {puts "Neither dnaA or polB positions could be identified for $species. Exiting"; exit 1}

			lappend ori_rep_tbl [join [concat $species $polB_pos] \t]
		} else {
			lappend ori_rep_tbl [join [concat $species $dnaA_pos] \t]
		}
	}
	return $ori_rep_tbl
}

# Get a list of the fully assembled Geobacillus genomes #
proc GetGeoAssembledSpecies {} {

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
	return $assembled_species_l
}

# Get the relative (to the origin) gene position of any gene in any assembled genomes #
proc GetRelativePosition {start_pos species_start relevant_genome_length} {
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
		set relative_pos [expr $start_pos - $ori_start]


		# If it's negative it lies downstream of the origin. Take the absolute value and leave as is #
		# If it's positive it lies upstream. We subtract it's position from the genome size. This will result in chromStart being > chromEnd #
		if {$relative_pos < 0} {
			set relative_pos [expr abs($relative_pos)]
		} else {
			set relative_pos [expr $relevant_genome_length - $relative_pos]
		}
	}
	
	## Get the fractional position ##
	set fraction_pos [tcl::mathfunc::roundto [expr double($relative_pos) / double($relevant_genome_length)] 8]
	return $fraction_pos
}

## 1. Set global variables ##
#############################

# List of reconciliation penalties at which genes in different sets (e.g. vertical/constant HGT/long HGT) will be considered #
set	penalty_list {3 4 5 6}
#set	penalty_list {6}
set	top_penalty [::tcl::mathfunc::max {*}$penalty_list]

# Data sets are the sets of genes for which we want to predict relative positions #
set data_sets {all vert const long}
#set data_sets {long}

# Open the sqlite ortholog database #
cd	/Users/aesin/Desktop/Geo_analysis/Geo_ortholog_nucl
sqlite3 db1 Geo_v_all_orth_location_NEW_db

# Log #
set		log_file /Users/aesin/desktop/Geo_analysis/HGT_position/Per_penalty_position.log
puts	stdout	"\nLog will be written to $log_file"
set 	out_log [open $log_file w]

## 2. Vital Geo paramaters ##
#############################

# Table of origin positions for each species #
multiputs stdout $out_log	"WORKING: Getting a table of the dnaA / dnaN (polB) gene positions to determine genome origin ..."
set		geo_ori_tbl [GetGeoOriPositions]
multiputs stdout $out_log	"DONE: Completed table of origin positions"

# List of fully assembled Geobacillus genomes #
multiputs stdout $out_log	"WORKING: Identifying fully assembled Geobacillus genomes ..."
set		assembled_geo_l [GetGeoAssembledSpecies]
multiputs stdout $out_log	"DONE: Identified fully assembled genomes"

# List of per species genome lengths (adjusted for plasmids) #
multiputs stdout $out_log	"WORKING: Getting a list of Geobacillus genomes lengths ..."
set		genome_length_l [split [string trim [openfile  /Users/aesin/Desktop/Geo_analysis/Geo_omes/Genome_lengths.tsv]] \n]
multiputs stdout $out_log	"DONE: Table of genome lengths"

## 3. Per penalty position data ##
##################################
# We use the penalty tip tables to see whether there are multiple transfer events per gene family in the consistent dataset
set per_penalty_tip_dir		/Users/aesin/Desktop/Mowgli/Mowgli_outputs/Per_penalty_tips
# The tip key translation tables are used to convert the per-group mowgli tip names back into unique protein IDs
set tip_key_dir				/Users/aesin/Desktop/Mowgli/Tip_keys_full

foreach penalty $penalty_list {

	# Receptor tables contain the tips for each transfer event in each gene family, along with the donor and receptor nodes #
	set penalty_tip_tbl_file $per_penalty_tip_dir/Per_penalty_tips_t$penalty\.tsv
	# File is new line separated and contains a header we remove #
	set penalty_tip_tbl [lrange [split [string trim [openfile $penalty_tip_tbl_file]] \n] 1 end]

	foreach set $data_sets {
		
		# Each penalty / set combination has a unique output path. Make a directory for each of the output types #
		set output_dir /users/aesin/desktop/Geo_analysis/HGT_position/For_circular/T$penalty/Per_species_$set
		file mkdir $output_dir/Full_entries

		# Counters and lists to keep track of getting gene positions #
		set event_counter 1

		set map_try_counter 0
		set map_fail_counter 0
		set plasmid_gene_counter 0

		set multi_xfer_const 0
		set total_geo_ids 0

		set plasmid_entries {}

		# For all genes, get all gene families for which we have at least one mapped Geobacillus gene #
		if {$set eq "all"} {
			set events [db1 eval {select DISTINCT group_number from t1 where binomial glob "Geobacillus*"}]
		} elseif {$set eq "vert"}	{
			if {$penalty == $top_penalty} {
				set vertical_file_name "Ver_const_t$penalty\.tsv"
			} else {
				set vertical_file_name "Ver_const_t$penalty\*$top_penalty\.tsv"	
			}
			set events [split [string trim [openfile [glob /Users/aesin/Desktop/Mowgli/Constant_events/Scenarios_1_2/$vertical_file_name]]] \n]
		} elseif {$set eq "const"}	{
			set const_HGT_file_name "HGT_const_t*$penalty\.tsv"
			set events [split [string trim [openfile [glob /Users/aesin/Desktop/Mowgli/Constant_events/Scenarios_1_2/$const_HGT_file_name]]] \n]
		} elseif {$set eq "long"}	{
			set long_dist_file_name "T$penalty\_full_lHGT_events.tsv"
			# First line is header #
			set events [lrange [split [string trim [openfile /Users/aesin/Desktop/Mowgli/Refined_events/Scenarios_1_2/Events/$long_dist_file_name]] \n] 1 end]
		}

		# For each assembled Geobacillus genome create a set of output tables
		foreach species $assembled_geo_l {
			set $species\_tbl {}
			set $species\_start_tbl {}
			set $species\_rel_start_tbl {}
		}

		set event_counter 1
		foreach event $events {

			# If this is a vertical event or any gene:
			if {$set eq "vert" || $set eq "all"} {
				# In these data, the event actually corresponds to unique groups
				set group $event

				# Get protein data from DB. Take only Geobacillus
				set entries_to_be_processed [db1 eval {SELECT prot_id from t1 where group_number = $group AND binomial glob "Geobacillus*";}]
				set total_geo_ids [expr $total_geo_ids + [llength $entries_to_be_processed]]
			}

			## If this is a transfer, we need to indentify a single unique event ##
			if {$set eq "long" || $set eq "const"} {

				# For lHGT input data structure see example of input file
				if {$set eq "long"} {
					set event_data_l [split $event \t]

					set group [lindex $event_data_l 0]
					set donor_edge [lindex $event_data_l 1]
					set receptor_edge [lindex $event_data_l 2]
					# The entries in this case are the mowgli tip names
					set entries_to_be_processed [split [lindex [split $event \t] 3] " "]

				# If const, then the event is just the group number 
				} else {
					set group $event
				}

				## Open the tip key translator file ##
				set tip_transl_tbl [split [string trim [openfile $tip_key_dir/TIP_KEY_$group\.tsv]] \n]

				

				if {$set eq "const"} {
					# Locate the corresponding entry in our table of receptor genomes
					set entry_in_recpt_tbl [lsearch -all -inline -glob $penalty_tip_tbl $group\t*]
					
					# If there is more than one HGT in this group, ignore it
					if {[llength $entry_in_recpt_tbl] > 1} {
						incr event_counter
						incr multi_xfer_const
						continue
					} elseif {[string length $entry_in_recpt_tbl] == 0} {
						puts "E1: not found $group"
						break 
					} else {
						set entries_to_be_processed [split [lindex [split [lindex $entry_in_recpt_tbl 0] \t] end] " "]
					}
				}
			}

			# `Entries` correspond to either protein IDs (from the all / vert sets) or mowgli tips involved in transfer (const / long sets)
			foreach entry $entries_to_be_processed {
				
				if {$set eq "long" || $set eq "const"} {
					set species_id_row [lsearch -all -inline -glob $tip_transl_tbl *\t$entry]

					if {[llength $species_id_row] != 1} {
						puts "Issue with finding the tip translation for group: $group and species tip: $entry"
						exit
					}

					set species_id_row_split [lindex [split [lindex $species_id_row 0] \t] 0] 
					set prot_id [string range $species_id_row_split [string last "\{" $species_id_row_split]+1 end-1]

				} elseif {$set eq "vert" || $set eq "all"} {
					set prot_id $entry
				}
				
				incr map_try_counter

				# Find the corresponding entry for the unique protein ID in the Geo_ortholog database #
				set prot_data_l [split [db1 eval {SELECT * from t1 where prot_id = $prot_id}] " "]
				if {[llength $prot_data_l] == 0} {
					incr map_fail_counter
					continue
				}

				## Remove any Anoxybacillus species and unassembled
				set species [lindex $prot_data_l 9]

				# Only continue if the genome is assembled #
				if {[lsearch $assembled_geo_l $species] == -1} {
					continue
				}

				################################################################################
				## Check the protein ID against the plasmid TAG list for the relevent species ##
				################################################################################

				## Open Plasmid tag file (if one exists for this species) ##
				set plasmid_tag_file /Users/aesin/Desktop/Geo_analysis/Geo_omes/Plasmid_tags/$species\_plasmid_locus_tags.tsv
				if {[file exists $plasmid_tag_file] == 1} {
					set tag_list [split [string trim [openfile $plasmid_tag_file]] \n]
					set transfer_locus [lindex $prot_data_l 8]

					## If this gene DOES map to the plasmid, do NOT add it to the total output list ##
					if {[lsearch $tag_list $transfer_locus] != -1} {
						incr plasmid_gene_counter
						lappend plasmid_entries $prot_data_l
						continue
					}
				}

				######################################################################################
				## Put the information (either full or just start position) into per-species tables ##
				######################################################################################

				## The strand will determine whether the true gene start position is the start or the end position as in the table. I.e. for complement strand, the "start" in the table is in fact the end, and vice versa ##
				set start_pos	[lindex $prot_data_l 5]
				set end_pos		[lindex $prot_data_l 6]

				###########################################################################
				## Make the distance to the origin a relative value of the genome length ##
				###########################################################################

				## Relative start positions ##
				if {[regexp "Geobacillus" $species] == 1} {
					set relevant_genome_length [lindex [split [lsearch -inline -glob $genome_length_l $species\t*] \t] 1]

					## Get the origin position for this genome ##
					set species_start [lsearch -glob -inline $geo_ori_tbl $species\t*]

					## If the origin lies on the complement strand, we take the "end" as the start. We also change the upstream/downstream orientation. ##
					set fraction_start_pos	[GetRelativePosition $start_pos $species_start $relevant_genome_length]
					set fraction_end_pos	[GetRelativePosition $end_pos $species_start $relevant_genome_length]

					# If the genome is not oriented relative to the origin, in some cases gene starts will be > gene ends. We reverse this #
					set ori_strand [lindex [split $species_start \t] 3]
					if {$ori_strand eq "-" && $fraction_end_pos < $fraction_start_pos} {
						set temp_start $fraction_start_pos
						set fraction_start_pos $fraction_end_pos
						set fraction_end_pos $temp_start
					}

					## Only makes sense to add branch data to the long dataset ##
					if {$set eq "long"} {
						## Append the fractional position to the full entry table ##
						lappend $species\_tbl [join [list [join $prot_data_l \t] [join $receptor_edge " "] [join $donor_edge " "] $fraction_start_pos $fraction_end_pos] \t]
						## Add fractional position to the per-species table ##
						lappend $species\_rel_start_tbl [join [list $group [join $receptor_edge " "] [join $donor_edge " "] $fraction_start_pos] \t]
					} else {
						## Append the fractional position to the full entry table ##
						lappend $species\_tbl [join [concat [join $prot_data_l \t] $fraction_start_pos $fraction_end_pos] \t]
						## Add fractional position to the per-species table ##
						lappend $species\_rel_start_tbl [join [concat $group $fraction_start_pos] " "]
					}
				} else {
					puts "Species is not Geobacillus. A1"
				}

			}

			puts "Penalty: $penalty == Set: $set == Events sorted: $event_counter / [llength $events]"
			incr event_counter

		}

		set combined_rel_pos {}
		foreach species $assembled_geo_l {
			set full_out_tbl [join [expr $$species\_tbl] \n]
			#set start_out_tbl [join [expr $$species\_start_tbl] \n]
			
			if {$set eq "long" || $set eq "const"} {
				set output_name "$species\_HGT.tsv"
			} else {
				set output_name "$species\_$set\.tsv"
			}

			set out [open $output_dir/Full_entries/$output_name w]
			puts $out $full_out_tbl
			close $out
			
			set rel_start_out_tbl [join [expr $$species\_rel_start_tbl] \n]
			lappend combined_rel_pos $rel_start_out_tbl
		}

		set out [open $output_dir/combined_rel_pos_assembled.tsv w]
		puts $out [join $combined_rel_pos \n]
		close $out
	 
	 	puts $out_log "\n\#\#\#\#\#\# Penalty: $penalty\tSet: $set \#\#\#\#\#\#"
		puts $out_log "Events processed: [expr $event_counter - 1]"
		if {$set eq "const"} {
			puts $out_log "Number of consistent transfer groups with > 1 transfer: $multi_xfer_const"
		} elseif {$set eq "vert" || $set eq "all"} {
			puts $out_log "Total geo ids processed: $total_geo_ids"
		}
		puts $out_log "Genes we tried to map: $map_try_counter"
		puts $out_log "Genes that failed to map: $map_fail_counter"
		puts $out_log "Genes mapping to plasmid: $plasmid_gene_counter"

	}
	puts $out_log "\n----------------------------"
}

close $out_log
db1 close

