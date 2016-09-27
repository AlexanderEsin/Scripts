#!/usr/local/bin/tclsh
## The script will isolate all the trees / groups that have HGT into Anoxy/Geobacillus coming from outside of the IG. If a group contains both transfers from OG and IG, it is conservatively not included. This is done at the most stringent penalty and then refined by checking whether the transfer still comes from an external source at all the more permissive thresholds ##

###########################################################################
## Procs ##
source ~/Dropbox/Scripts/General_utils.tcl

## Paths ##
set direct /users/aesin/desktop/Mowgli
set mow_outs_direct $direct/Mowgli_outputs

set input_direct $direct/Consistent_HGT_Vertical
set output_direct $direct/Long_distance_HGT

set inside_included TRUE
set all_scenarios FALSE

if {$inside_included == TRUE} {
	if {$all_scenarios == TRUE} {
		set input_path  $input_direct/Consistent_Full/All_scenarios
		set output_direct $output_direct/Full/All_scenarios
	} else {
		set input_path  $input_direct/Consistent_Full/Scenarios_1_2
		set output_direct $output_direct/Full/Scenarios_1_2
	}
} else {
	if {$all_scenarios == TRUE} {
		set input_path  $input_direct/Consistent_NO_IG/All_scenarios
		set output_direct $output_direct/NO_IG/All_scenarios
	} else {
		set input_path  $input_direct/Consistent_NO_IG/Scenarios_1_2
		set output_direct $output_direct/NO_IG/Scenarios_1_2
	}
}

set output_direct $output_direct/Per_penalty_refined
file mkdir $output_direct

## List of penalties for which to perform long-distance transfer refinement ##
set refined_transfer_penalty_list [list 4 5 6 8 10 20]

###########################################################################

## Get the list of inside group nodes (result of Internal_nodes_IG.R) ##
set internal_IG_nodes_l [split [string trim [openfile $direct/Inside_group_nodes.txt]] \n]

######################################

foreach orig_penalty $refined_transfer_penalty_list {

	puts "\n######################################\nTesting for penalty $orig_penalty :"

	## Set the penalties at which we will see whether the transfer is external ##
	set full_penalty_list	[list 3 4 5 6 8 10 20]
	set tested_penalty_list	[lrange $full_penalty_list 0 [lsearch $full_penalty_list $orig_penalty]]
	set reversed_test_list	[reverse_dict $tested_penalty_list]

	######################################
	## Open the list of HGTs at the top penalty tested ##
	set penalty_list_str "t[join $tested_penalty_list "_t"]"
	set HGT_tree_l [split [string trim [openfile $input_path/true_HGT_intersect_$penalty_list_str\.tsv]] \n]

	## A list holder for those groups that have all HGTs from an external source across all penalties ##
	set counter 0

	######################################

	foreach penalty $reversed_test_list {
		puts "\n\tAt penalty: $penalty"
		set total_HGT_events 0
		set from_IG 0
		set external_HGT_only 0

		## At every iteration apart from the first, use the refined group list from the previous iteration at input for testing ##
		if {$counter != 0} {
			set HGT_tree_l $external_list
		}
		set external_list {}


		foreach HGT_candidate $HGT_tree_l {
			set from_external 1

			## Open the transfers file to get a list of all the HGT events into Anoxy/Geobacillus for this group ##
			set HGT_events [regexp -all -line -inline {(?:HGT)+.+} [openfile $mow_outs_direct/Mow_test_out_t$penalty/$HGT_candidate/Transfers_in.tsv]]
			set total_HGT_events [expr $total_HGT_events + [llength $HGT_events]]

			## For each event, identify the donor 
			foreach event $HGT_events {
				set donor_edge [string range $event [string first " " $event]+1 [string first \t $event]-1]
				set donor_nodes [split $donor_edge \,]
				if {[lsearch $internal_IG_nodes_l [lindex $donor_nodes 0]] != -1 || [lsearch $internal_IG_nodes_l [lindex $donor_nodes 1]] != -1} {
					# puts "From internal: $HGT_candidate"
					incr from_IG
					set from_external 0
				}
			}
			if {$from_external == 1} {
				incr external_HGT_only
				lappend external_list $HGT_candidate
			}
		}

		incr counter

		puts "\tTotal trees: [llength $HGT_tree_l]"
		puts "\tTotal HGT events: $total_HGT_events"
		puts "\tTotal events from IG donor: $from_IG"
		puts "\tTrees with all HGT from OUTSIDE IG: $external_HGT_only"
	}

	set out [open $output_direct/T$orig_penalty\_refined_groups.tsv w]
	puts $out [join $external_list \n]
	close $out
}






