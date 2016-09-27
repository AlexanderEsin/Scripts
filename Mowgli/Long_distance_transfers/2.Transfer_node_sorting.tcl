#!/usr/local/bin/tclsh
source ~/Dropbox/Scripts/General_utils.tcl
package require sqlite3

###########################################################################
set direct /users/aesin/desktop/Mowgli

## Pick penalty set and minimum sequences in group for count to be included ##
set test_penalty_list 	[list 4 5 6 8 10 20]
set test_min_taxa_list	[list 0 50]

## Other options ##
set inside_included TRUE
set all_scenarios FALSE

## Set up directories
set input_dir $direct/Long_distance_HGT
set db_dir /users/aesin/desktop/Geo_analysis/Geo_v_all/2.0/Geo_ortholog_db

if {$inside_included == TRUE} {
	if {$all_scenarios == TRUE} {
		set input_dir  $input_dir/Full/All_scenarios
	} else {
		set input_dir  $input_dir/Full/Scenarios_1_2
	}
} else {
	if {$all_scenarios == TRUE} {
		set input_dir  $input_dir/NO_IG/All_scenarios
	} else {
		set input_dir  $input_dir/NO_IG/Scenarios_1_2
	}
}

set output_dir $input_dir
set input_dir $input_dir/Per_penalty_refined

###########################################################################
## Make some output folders ##
file mkdir $output_dir/Receptor_edges
file mkdir $output_dir/Donor_edges
file mkdir $output_dir/At_root/AG_root
file mkdir $output_dir/At_root/G_root

###########################################################################

foreach min_taxa $test_min_taxa_list {

	foreach penalty $test_penalty_list {

		puts "\n######################################\nTesting for penalty $penalty and min. sequences included $min_taxa :"

		## If we need to filter by number of sequences, open the relevant database ##
		if {$min_taxa > 4} {
			sqlite3 db1 $db_dir
			set output_string_addon "_min_$min_taxa"
		} else {
			set output_string_addon ""
		}

		#####################################
		## Open the list of refined groups ##
		set refined_group_l [split [string trim [openfile $input_dir/T$penalty\_refined_groups.tsv]] \n]
		set total_refined_groups [llength $refined_group_l]

		set total_HGT_events 0
		set total_groups_above_min_taxa 0

		set receptor_edge_l {}
		set donor_edge_l {}
		set receptor_donor_l {}

		set geobac_root_groups_l {}
		set at_root_l {}

		#####################################
		foreach group $refined_group_l {

			## If we have a minimum taxa threshold, check how many sequences we have in the group, and ignore any groups that have fewer sequences ##
			if {$min_taxa > 4} {
				set taxa [db1 eval {select count(*) from t1 where group_num = $group}]
				if {$taxa < $min_taxa} {
					continue
				} else {
					incr total_groups_above_min_taxa
				}
			}

			## Pull out all HGT events for the group for the penalty threshold. We are alreay certain that all the transfers are from outside the IG for this penalty and all more permissive penalties (see External_transfers_only.tcl) ##
			set HGT_events [regexp -all -line -inline {(?:HGT)+.+} [openfile $direct/Mowgli_outputs/Mow_test_out_t$penalty/$group/Transfers_in.tsv]]
			set total_HGT_events [expr $total_HGT_events + [llength $HGT_events]]
			
			# puts $HGT_events

			foreach event $HGT_events {
				set receptor_edge [string range $event [string last \t $event]+1 end]
				set receptor_nodes [split $receptor_edge \,]

				set donor_edge [string trim [string range $event [string first \: $event]+1 [string last \t $event]]]
				set donor_nodes [split $donor_edge \,]
				## Let's add the group number to the donor edges so that we can refer to the taxa present in that group specifically ##
				set donor_nodes [concat $group $donor_nodes]

				if {$receptor_edge eq "3643,3642"} {
					lappend geobac_root_groups_l $group
				} elseif {[lindex $receptor_nodes 0] eq "3728"} {
					lappend at_root_l $group
				}

				lappend receptor_edge_l [join $receptor_nodes \t]
				lappend donor_edge_l [join $donor_nodes \t]
				lappend receptor_donor_l "[join $receptor_nodes " "]\t[join $donor_nodes " "]"
			}
		}

		puts "\n\tTotal groups: $total_refined_groups"
		if {$min_taxa > 4} {puts "\t\tTotal groups with $min_taxa sequences or more: $total_groups_above_min_taxa"; db1 close}
		puts "\tTotal HGT events: $total_HGT_events"
		puts "\tReceptor edges found: [llength $receptor_edge_l]"
		puts "\tDonor edges found: [llength $donor_edge_l]"

		# puts [llength $geobac_root_groups_l]

		#####################################
		## Output all receptor edges ##
		set out [open $output_dir/Receptor_edges/T$penalty\_refined_receptor_edges$output_string_addon\.txt w]
		puts $out [join $receptor_edge_l \n]
		close $out

		## Output all donor edges ##
		set out [open $output_dir/Donor_edges/T$penalty\_refined_donor_edges$output_string_addon\.txt w]
		puts $out [join $donor_edge_l \n]
		close $out

		## Output the donor edges with the receptor edges. Format: Receptor_edge\tDonor_edge ##
		set out [open $output_dir/Donor_edges/T$penalty\_refined_receptor_donor_list$output_string_addon\.txt w]
		puts $out [join $receptor_donor_l \n]
		close $out


		## Output list of groups where the receptor edge is the Geobacillus root ##
		set out [open $output_dir/At_root/G_root/T$penalty\_refined_receptor_edges$output_string_addon\_geobac_root.txt w]
		puts $out [join $geobac_root_groups_l \n]
		close $out

		## Output list of groups where the receptor edge is the Anoxy/Geobacillus root ##
		set out [open $output_dir/At_root/AG_root/T$penalty\_refined_receptor_edges$output_string_addon\_all_root.txt w]
		puts $out [join $at_root_l \n]
		close $out

		#####################################

	}
}




