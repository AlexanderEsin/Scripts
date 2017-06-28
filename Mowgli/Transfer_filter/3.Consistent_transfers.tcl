#!/usr/local/bin/tclsh
## The script will isolate all the trees / groups that have HGT into Anoxy/Geobacillus coming from outside of the IG. For each gene family (group) we get a list of all the transfers into Anoxy/Geobacillus. If the transfer is predicted to come from a long-distance donor (i.e. not IG), we also check whether the Anoxy/Geo species involved at that penalty (e.g. 6) are also the receptors of a long distance transfer in more permissive penalties (i.e. 5, 4, 3). Thus, we only identify consistent long-distance transfers in each gene family. In some cases, a single gene family will contain transfers that are from IG donors and some from long-distance donors; this script differentiates the two, and so will include gene families which still containing consistent long-distance HGTs even if there are HGT from IG groups also involved. Thus, the output is a list of groups and associated tips which show consistent long-distance HGT across a range of penalties. ##

## 0. Procs ##
##############
source ~/Documents/Scripts/General_utils.tcl

proc GetDonorNodes {group_hgt_event} {
	set columns		[split [string trim $group_hgt_event] \t]
	if {[llength $columns] != 4} {puts "Why not 4 columns in this event?: $group_hgt_event \{GetDonorNodes\}"; break 1}
	set donor_edge	[lindex [split [lindex $columns 1] " "] 1]
	set donor_nodes	[split $donor_edge ","]
	return $donor_nodes
}

proc GetReceptorNodes {group_hgt_event} {
	set columns		[split [string trim $group_hgt_event] \t]
	if {[llength $columns] != 4} {puts "Why not 4 columns in this event?: $group_hgt_event \{GetReceptorNodes\}"; break 1}
	set receptor_edge	[lindex [split [lindex $columns 1] " "] end]
	set receptor_nodes	[split $receptor_edge ","]
	return $receptor_nodes
}

proc GetHgtTips {group_hgt_event} {
	set columns		[split [string trim $group_hgt_event] \t]
	if {[llength $columns] != 4} {puts "Why not 4 columns in this event? $group_hgt_event \{GetHgtTips\}"; break 1}
	set HGT_tips [split [lindex $columns end] " "]
	return $HGT_tips
}

proc GetTransferType {donor_nodes IG_nodes} {
	if {[llength $donor_nodes] != 2} {puts "There must be two donor nodes (one edge = 2 nodes)"; exit 1}

	if {[lsearch $IG_nodes [lindex $donor_nodes 0]] != -1 || [lsearch $IG_nodes [lindex $donor_nodes 1]] != -1} {
		set transfer_type internal
	} else {
		set transfer_type external
	}
	return $transfer_type
}

proc GetConsistentEvents {events} {
	
	# Define global variables in two calls - first is normal variables. Second is a double dereference referring to events at a specific penalty
	global test_penalty HGT_group IG_nodes_l per_group_log per_group_hgt_loss per_group_hgt_inconsistent
	global $test_penalty\_tip_data

	multiputs $per_group_log stdout "\nPenalty: $test_penalty"

	# Get the set of events at a given penalty, against which we will be testing the transfer events
	set test_events [lsearch -glob -inline -all [subst $$test_penalty\_tip_data] $HGT_group\t*]

	# For each transfer event in the group...
	foreach event $events {
		# Get event and the type (internal or external)
		set HGT_event	[string trim [lindex [split $event {|}] 0]]
		set HGT_type	[string trim [lindex [split $event {|}] 1]]

		# Get the event tips to be tested
		set HGT_event_tips [GetHgtTips $HGT_event]
		multiputs $per_group_log stdout "\n$HGT_type event: $HGT_event"

		foreach HGT_tip $HGT_event_tips {
			set HGT_event_contains_tip [lsearch -glob -inline -all $test_events *$HGT_tip*]
			# It's possibe that the events at a lower penalty no longer contain the tips from events at higher penalties. I.e. HGT is no longer predicted into those tips - however, there may still be a predicted HGT into Geobacillus elsewhere in the gene tree
			if {[llength $HGT_event_contains_tip] == 0} {
				set events [lremove $events $event]

				multiputs $per_group_log stdout "Tip no longer part of an HGT. Removing"
				lappend per_group_hgt_loss $HGT_event
				break
			}

			set next_level_HGT_type [GetTransferType [GetDonorNodes $HGT_event_contains_tip] $IG_nodes_l]
			multiputs $per_group_log stdout "$HGT_tip. Transfer type == $next_level_HGT_type"

			if {$next_level_HGT_type != $HGT_type} {
				set events [lremove $events $event]

				multiputs $per_group_log stdout "Tip $HGT_tip at penalty $test_penalty is not an $HGT_type event. Adding to inconsistent events..."
				lappend per_group_hgt_inconsistent $HGT_event
				break
			} else {
				multiputs $per_group_log stdout "Tip $HGT_tip is still part of $HGT_type event at $test_penalty"
			}
		}
	}

	return $events
}

## 1. Global variables ##
#########################

# Paths and output folders
set mowgli_dir /Users/aesin/Desktop/Mowgli

# List of penalties for which to perform long-distance transfer refinement
set long_HGT_penalty_l [list 3 4 5 6]

## 2. I/O folder assignment ##
##############################

# Scenarios BS - need to revisit. I/O paths
set all_scenarios FALSE

if {$all_scenarios == TRUE} {
	set input_path  $mowgli_dir/Constant_events/All_scenarios
	set output_path $mowgli_dir/Refined_events/All_scenarios
} else {
	set input_path  $mowgli_dir/Constant_events/Scenarios_1_2
	set output_path $mowgli_dir/Refined_events/Scenarios_1_2
}

# Make the necessary ouput folders #
foreach penalty $long_HGT_penalty_l {
	file mkdir $output_path/Per_group/$penalty
}
file mkdir $output_path/Groups
file mkdir $output_path/Events


## 3. Read in extra data ##
###########################

# Read in the per-penalty group/node/tip data #
set per_penalty_HGT_folder $mowgli_dir/Mowgli_outputs/Per_penalty_tips
foreach penalty $long_HGT_penalty_l {
	set $penalty\_tip_data [lrange [split [string trim [openfile $per_penalty_HGT_folder/Per_penalty_tips_t$penalty\.tsv]] \n] 1 end]
}

# Get the list of inside group nodes (result of 0. Internal_nodes_IG.R)
set IG_nodes_l [split [string trim [openfile $mowgli_dir/Inside_group/Inside_group_nodes.txt]] \n]

## 4. Per penalty lHGT prediction ##
####################################

set cross_pen_log [open $output_path/Full_log.log w]

foreach orig_penalty $long_HGT_penalty_l {
	multiputs $cross_pen_log stdout "\n######################################\nTesting for penalty $orig_penalty :"

	# We test top down 
	set test_penalty_l	[lrange $long_HGT_penalty_l 0 [lsearch $long_HGT_penalty_l $orig_penalty]]
	set reverse_pen_l	[lrange [reverse_dict $test_penalty_l] 1 end]

	multiputs $cross_pen_log stdout "Long distant transfer consistent over penalties: [join $reverse_pen_l " "]"
	
	# List of consistent HGT groups at the penalty of interest #
	set penalty_list_str "t[join $test_penalty_l "_t"]"
	set HGT_groups_l [split [string trim [openfile $input_path/HGT_const_$penalty_list_str\.tsv]] \n]
	
	set constHGT_group_l {}
	set lHGT_group_l {}
	set sHGT_group_l {}

	set per_penalty_lHGT_out [list "Group\tDonor_nodes\tReceptor_nodes\tTips"]
	set per_penalty_sHGT_out [list "Group\tDonor_nodes\tReceptor_nodes\tTips"]

	set per_penalty_inconsistent {}
	set per_penalty_loss {}

	foreach HGT_group $HGT_groups_l {

		set per_group_log [open $output_path/Per_group/$orig_penalty/$HGT_group\.txt w]
		multiputs $per_group_log stdout "#####\nGroup: $HGT_group"

		# Total HGT events for this group at this penalty:
		set independent_events [lsearch -glob -inline -all [expr $$orig_penalty\_tip_data] $HGT_group\t*]
		

		set group_hgt_events {}
		
		set per_group_hgt_loss {}
		set per_group_hgt_inconsistent {}

		foreach event $independent_events {
			set donor_nodes [GetDonorNodes $event]
			set transfer_type [GetTransferType $donor_nodes $IG_nodes_l]

			lappend group_hgt_events [join [list $event $transfer_type] "\t|\t"]
			multiputs $per_group_log stdout "$event. Transfer type == $transfer_type"
		}

		## If we are testing for consistency ..
		if {[llength $reverse_pen_l] > 0} {
			foreach test_penalty $reverse_pen_l {

				if {[llength $group_hgt_events] == 0} {
					multiputs $per_group_log stdout "\nNo tips came from a consistent HGT source in the above penalty."
					break
				}
				puts $group_hgt_events
				set group_hgt_events [GetConsistentEvents $group_hgt_events]				
			}

			if {[llength $group_hgt_events] > 0} {
				lappend constHGT_group_l $HGT_group

				foreach constEvent $group_hgt_events {
					set HGT_event	[string trim [lindex [split $constEvent {|}] 0]]
					set HGT_type	[string trim [lindex [split $constEvent {|}] 1]]

					# Format the output data
					set donor_nodes [join [GetDonorNodes $HGT_event] " "]
					set receptor_nodes [join [GetReceptorNodes $HGT_event] " "]
					set tips [join [GetHgtTips $HGT_event] " "]
					set out_line_str [join [list $HGT_group $donor_nodes $receptor_nodes $tips] \t]

					# Append the output to the relevant short or long transfer list
					if {$HGT_type eq "internal"} {
						lappend per_penalty_sHGT_out $out_line_str
						lappend sHGT_group_l $HGT_group
					} else {
						lappend per_penalty_lHGT_out $out_line_str
						lappend lHGT_group_l $HGT_group
					}
				}
			} else {
				multiputs $per_group_log stdout "\nNo consistent short or long distance HGT events in $HGT_group"
			}

			multiputs $per_group_log stdout "\nEvents that are inconsistently long/short transfers:"
			if {[llength $per_group_hgt_inconsistent] != 0} {
				multiputs $per_group_log stdout [join $per_group_hgt_inconsistent \n]
				lappend per_penalty_inconsistent $per_group_hgt_inconsistent
			}
			
			multiputs $per_group_log stdout "\nHGT events that are no longer HGTs:"
			if {[llength $per_group_hgt_loss] != 0} {
				multiputs $per_group_log stdout [join $per_group_hgt_loss \n]
				lappend per_penalty_loss $per_group_hgt_loss
			}

		## Otherwise, we just write out the results for the individual penalty ##
		} else {
			foreach hgt_event $group_hgt_events {
				set HGT_event	[string trim [lindex [split $hgt_event {|}] 0]]
				set HGT_type	[string trim [lindex [split $hgt_event {|}] 1]]

				# Format the output data
				set donor_nodes [join [GetDonorNodes $HGT_event] " "]
				set receptor_nodes [join [GetReceptorNodes $HGT_event] " "]
				set tips [join [GetHgtTips $HGT_event] " "]
				set out_line_str [join [list $HGT_group $donor_nodes $receptor_nodes $tips] \t]

				# Append the output to the relevant short or long transfer list
				if {$HGT_type eq "internal"} {
					lappend per_penalty_sHGT_out $out_line_str
					lappend sHGT_group_l $HGT_group
				} else {
					lappend per_penalty_lHGT_out $out_line_str
					lappend lHGT_group_l $HGT_group
				}
			}
		}

		close $per_group_log
		
	}

	set unique_lHGT_groups [lsort -unique $lHGT_group_l]
	set unique_sHGT_groups [lsort -unique $sHGT_group_l]


	# Write out full events
	set out [open $output_path/Events/T$orig_penalty\_full_lHGT_events.tsv w]
	puts $out [join $per_penalty_lHGT_out \n]
	close $out

	set out [open $output_path/Events/T$orig_penalty\_full_sHGT_events.tsv w]
	puts $out [join $per_penalty_sHGT_out \n]
	close $out

	# Write out just the groups
	set out [open $output_path/Groups/T$orig_penalty\_unique_lHGT_groups.tsv w]
	puts $out [join $unique_lHGT_groups \n]
	close $out
	
	set out [open $output_path/Groups/T$orig_penalty\_unique_sHGT_groups.tsv w]
	puts $out [join $unique_sHGT_groups \n]
	close $out

	multiputs $cross_pen_log stdout		"\nOriginal num of HGT groups: [llength $HGT_groups_l]\
										\nGroups with consistent long and/or short HGTS: [llength $constHGT_group_l]\
										\n\tGroups with consistent long-distance HGT tips: [llength $unique_lHGT_groups]\
										\n\tGroups with consistent short-distance HGT tips: [llength $unique_sHGT_groups]\
										\n\nIndependent lHGT events: [expr [llength $per_penalty_lHGT_out] -1]\
										\nIndependent sHGT events: [expr [llength $per_penalty_sHGT_out] -1]"
}

close $cross_pen_log
puts "\n\nDONE\n"
puts "Output and log written to: $output_path"




