proc GetDonorNodes {group_hgt_event} {
	set columns		[split [string trim $group_hgt_event] \t]
	if {[llength $columns] != 6} {puts "Why not 6 columns in this event?: $group_hgt_event \{GetDonorNodes\}"; exit 1}
	# Get the edge and nodes
	set donor_edge	[string range [lindex $columns 1] 1 end-1]
	set donor_nodes	[split $donor_edge ","]
	return $donor_nodes
}

proc GetReceptorNodes {group_hgt_event} {
	set columns		[split [string trim $group_hgt_event] \t]
	if {[llength $columns] != 6} {puts "Why not 4 columns in this event?: $group_hgt_event \{GetReceptorNodes\}"; exit 1}
	# Get the edge and nodes
	set receptor_edge	[string range [lindex $columns 2] 1 end-1]
	set receptor_nodes	[split $receptor_edge ","]
	return $receptor_nodes
}

proc GetHgtTips {group_hgt_event} {
	set columns		[split [string trim $group_hgt_event] \t]
	if {[llength $columns] != 6} {puts "Why not 4 columns in this event? $group_hgt_event \{GetHgtTips\}"; exit 1}
	set HGT_tips [split [string trim [lindex $columns end]] " "]
	return $HGT_tips
}

proc GetTransferType {donor_nodes IG_nodes} {
	if {[llength $donor_nodes] != 2} {puts "There must be two donor nodes (one edge = 2 nodes)"; exit  1}

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
		puts "\n$event"
		# Get event and the type (internal or external)
		set HGT_event	[string trim [lindex [split $event {|}] 0]]
		set HGT_type	[string trim [lindex [split $event {|}] 1]]

		# Get the event tips to be tested
		set HGT_event_tips [GetHgtTips $HGT_event]
		multiputs $per_group_log stdout "\n$HGT_type event: $HGT_event"

		foreach HGT_tip $HGT_event_tips {
			# puts "\n\n$HGT_tip"
			set HGT_event_contains_tip [lsearch -glob -inline -all $test_events "* $HGT_tip *"]
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
