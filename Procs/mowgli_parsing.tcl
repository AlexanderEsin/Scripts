## Procs for Mowgli output parsing ##

proc ParseMappingFile {mapping_file} {
	# Read in the Fullmapping file 
	set recon_map [string trim [openfile $mapping_file]]

	# If empty return empty events list
	if {[string length $recon_map] == 0} {
		set events_list {}
		return $events_list
	}

	# Remove extra tabs and newlines
	regsub -all {\n\n} $recon_map "\n" recon_map
	regsub -all {\t\t} $recon_map "\t" recon_map

	# Remove header (first two lines after whitespace trim)
	set events_list [lrange [split $recon_map \n] 2 end]
	return $events_list
}

proc GetTransferEvents {events} {
	set transfer_events {}
	foreach event $events {
		if {[lsearch -glob $event *Tran*] > 0} {
			lappend transfer_events $event
		}
	}
	return $transfer_events
}

proc GetFinalSubevent {event} {
	set event_data [split [string trim $event] \t]
	set final_subevent [lindex $event_data end]
	return $final_subevent
}

proc GetEventType {subevent} {
	set event_type [lindex [split [string trim $subevent] " "] end]
	return $event_type
}

proc IsSubeventTransfer {subevent} {
	set is_transfer FALSE
	if {[regexp {;} $subevent] == 1} {
		set is_transfer TRUE
	}
	return $is_transfer
}

proc IsSubeventLoss {subevent} {
	set is_loss FALSE
	if {[regexp {Loss} $subevent] == 1} {
		set is_loss TRUE
	}
	return $is_loss
}

proc PutLossAtBack {non_trans_events} {

	if {[llength $non_trans_events] > 1} {
		foreach event $non_trans_events {
			set final_subevent [GetFinalSubevent $event]
			if {[IsSubeventLoss $final_subevent] == TRUE} {
				set non_trans_events [lremove $non_trans_events $event]
				lappend non_trans_events $event
			}
		}
	}

	return $non_trans_events
}

proc ReduceRedundantEdges {redundant_events} {
	set non_redun_out {}

	foreach event $redundant_events {
		set final_subevent	[GetFinalSubevent $event]
		set parsed_event	[ParseNontransSubevent $final_subevent]
		set edge			[dict get $parsed_event edge]

		if {[lsearch -glob $non_redun_out *$edge*] == -1} {
			lappend non_redun_out $event
		}
	}
	return $non_redun_out
}

proc ReduceTransferBranch {receiver_edge event_type events_list dir_log} {
	set tested_edge $receiver_edge

	while 1 {

		# If the transfer is into a single taxon - there can be no lower node
		if {$event_type eq "Extant"} {
			break
		}

		set transfered_edge_child [lindex [split $tested_edge \,] 1]
		set children [lsearch -all -inline -glob $events_list *\($transfered_edge_child,*]

		set branch_fates {}
		foreach child_event $children {
			set final_subevent [GetFinalSubevent $child_event]
			if {[IsSubeventTransfer $final_subevent] == TRUE || [IsSubeventLoss $final_subevent] == TRUE} {
				continue
			} else {
				lappend branch_fates $final_subevent
			}
		}

		## If there are duplications, for each duplication picked up - there will be two further speciations - we can reduce this down to a single entry for the fate ##
		set branch_fates [ReduceRedundantEdges $branch_fates]
	
		## Check whether both speciated, or one was lost ##
		if {[llength $branch_fates] == 2} {
			multiputs $dir_log stdout "This is the final edge: $tested_edge"
			set receiver_edge $tested_edge
			break
		} elseif {[llength $branch_fates] == 1} {
			set parsed_event	[ParseNontransSubevent [lindex $branch_fates 0]]
			set event_type		[dict get $parsed_event event_type]
			set edge			[dict get $parsed_event edge]

			# If of two fates, one is loss and the other is extant, there can be no further down-sampling of the tree #
			if {$event_type eq "Extant"} {
				set receiver_edge [string range $edge 1 end-1]
				multiputs $dir_log stdout "This is the final edge: $receiver_edge"
				break
			} else {
				multiputs $dir_log stdout "One branch was lost, continuing down..."
				set tested_edge [string range $edge 1 end-1]
			}
		} else {
			multiputs $dir_log stdout "Error: loss of both fates or more than two fates for branch: Tree: $out_dir A2"
			exit
		}
	}

	return $receiver_edge
}

proc ParseEvent {event} {

	set split_event [split [string trim $event] \t]

	set genet_edge [lindex $split_event 0]
	set parsed_event [dict create genet_edge $genet_edge]

	# Parsing the gene tree nodes is agnostic to brackets because the root event uses (x,x) and all other events use {x,x}
	dict set parsed_event genet_parent	[string range [lindex [split $genet_edge \,] 0] 1 end]
	dict set parsed_event genet_child	[string range [lindex [split $genet_edge \,] 1] 0 end-1]

	set final_subevent [lindex $split_event end]
	if {[IsSubeventTransfer $final_subevent] == TRUE} {
		set parsed_subevent [ParseTransferSubevent $final_subevent]
		dict set parsed_event type "Trans"
	} else {
		set parsed_subevent [ParseNontransSubevent $final_subevent]
		dict set parsed_event type "NonTrans"
	}

	set parsed_event [dict merge $parsed_event $parsed_subevent]
	return $parsed_event
}

proc ParseTransferSubevent {subevent} {

	set subevent_split [split [string range $subevent 1 [string first " " $subevent]-2] \;]

	set transfer_parse [dict create event_type [string trim [string range $subevent [string first " " $subevent] end]]]

	dict set transfer_parse genet_xfer_child		[lindex $subevent_split 0]
	dict set transfer_parse donor_edge				[lindex $subevent_split 1]
	dict set transfer_parse receiver_edge			[lindex $subevent_split 2]
	
	return $transfer_parse
}

proc ParseNontransSubevent {subevent} {
	set subevent_split [split $subevent " "]
	set edge_split [split [string range [lindex $subevent_split 0] 1 end-1] \,]

	set event_parse [dict create event_type [lindex $subevent_split 1]]
	dict set event_parse edge	[lindex $subevent_split 0]
	dict set event_parse parent	[lindex $edge_split 0]
	dict set event_parse child	[lindex $edge_split 1]

	return $event_parse
}

proc GetCorrectChild {child_events key_edge} {
	foreach child_event $child_events {
		set parsed_child	[ParseEvent $child_event]
		if {[dict get $parsed_child type] eq "Trans"} {
			set tested_edge	[dict get $parsed_child donor_edge]
		} else {
			set tested_edge	[dict get $parsed_child edge]
		}

		if {$tested_edge != $key_edge} {
			set child_events [lremove $child_events $child_event]
		}
	}
	return $child_events
}

proc IsTransIntoEdge {event test_edge} {
	set TransIn FALSE

	set parsed_event	[ParseEvent $event]
	set receiver_edge	[dict get $parsed_event receiver_edge]

	if {$receiver_edge == $test_edge} {
		set TransIn TRUE
	}

	return $TransIn
}

proc SkipAllTransOut {trans_event test_edge events_list {direction down}} {
	set next_event $trans_event

	set parsed_event [ParseEvent $next_event]
	set type [dict get $parsed_event type]

	# This loop has to run at least once because the event provided must be a transfer OUT #
	while {$type eq "Trans" && [IsTransIntoEdge $next_event $test_edge] == FALSE} {
		set previous_event $next_event

		# For every transfer up - there can only be one parent
		if {$direction eq "up"} {
			set genet_parent	[dict get $parsed_event genet_parent]
			set next_event		[lsearch -all -inline -glob $events_list "\{*,$genet_parent\}*"]
			# We might skip up to the root
			if {[llength $next_event] == 0} {
				set next_event		[lsearch -all -inline -glob $events_list "\(*,$genet_parent\)*"]
			}
		# For every transfer down - there will be two children
		} elseif {$direction eq "down"} {
			set genet_child		[dict get $parsed_event genet_child]
			set child_events	[lsearch -all -inline -glob $events_list "\{$genet_child\,*"]		
			set next_event		[GetCorrectChild $child_events $test_edge]
		} else {
			error "Direction can only be \"up\" or \"down\""
			exit
		}
		
		if {[llength $next_event] != 1} {
			error "Expecting to find only one event. Code 3451"
			exit
		}
		set next_event		[lindex $next_event 0]
		set parsed_event	[ParseEvent $next_event]
		set type			[dict get $parsed_event type]
	}
	if {$direction eq "up"} {
		set skipped_events [dict create parent_event $next_event]
		dict set skipped_events penul_event $previous_event
	} else {
		set skipped_events [dict create child_event $next_event]
		dict set skipped_events penul_event $previous_event
	}
	
	return $skipped_events
}

proc GetParentEvent {child_event events_list} {
	set parsed_child [ParseEvent $child_event]
	set child_genet_parent	[dict get $parsed_child genet_parent]
	set child_edge			[dict get $parsed_child edge]

	set parent_event		[lsearch -all -inline -glob $events_list \{*,$child_genet_parent\}*]

	# Could be root - it has a different format: ( ,112)	(3640,3639) Spec0
	if {[llength $parent_event] == 0} {
		set parent_event		[lsearch -all -inline -glob $events_list "\( ,$child_genet_parent\)\t*"]
	}

	# If we still can't find a parent event - error.
	if {[llength $parent_event] == 0} {
		error "Expecting one parent event. Code: D0GS"
	}

	set parent_event		[lindex $parent_event 0]
	set parsed_parent		[ParseEvent $parent_event]

	if {[dict get $parsed_parent type] eq "Trans" && [IsTransIntoEdge $parent_event $child_edge] == FALSE} {
		set parent_child_events [SkipAllTransOut $parent_event $child_edge $events_list "up"] 
		dict set parent_child_events child_event $child_event
		puts "SKIPPED TRANSFER UP: [dict get $parent_child_events parent_event]"
	} else {
		set parent_child_events [dict create parent_event $parent_event]
		dict set parent_child_events penul_event $child_event
		dict set parent_child_events child_event $child_event
	}
	return $parent_child_events
}

proc GetSisterEvent {parent_event known_child events_list} {
	set parsed_parent		[ParseEvent $parent_event]
	set parent_genet_child	[dict get $parsed_parent genet_child]

	set child_events [lsearch -all -inline -glob $events_list \{$parent_genet_child,*\}\t*]
	if {[llength $child_events] != 2} {
		error "There should be two children??. Code 1190"
	}
	
	set sister_event [lremove $child_events $known_child]
	if {[llength $sister_event] != 1} {
		error "There should be one sister??. Code 1195"
	}

	return [lindex $sister_event 0]
}

proc TestAllChildren {start_event events_list} {
	set events_to_test [list $start_event]

	set condition "Loss"
	while {[llength $events_to_test] != 0} {

		set tested_event	[lindex $events_to_test 0]
		puts $tested_event

		set event_info		[ParseEvent $tested_event]
		set type			[dict get $event_info type]

		if {$type eq "Trans"} {
			set test_edge			[dict get $event_info donor_edge]
			set child_parent_events	[SkipAllTransOut $tested_event $test_edge $events_list "down"]
			set next_event			[dict get $child_parent_events child_event]
			set event_info			[ParseEvent $next_event]
			puts "SKIPPED TRANSFER DOWN: $next_event"
		}

		set event_type		[dict get $event_info event_type]

		# If the event is a loss, remove and continue. If extant, return success
		if {$event_type eq "Loss"} {
			set events_to_test [lremove $events_to_test $tested_event]
			continue
		} elseif {$event_type eq "Extant"} {
			set condition "Ext: $tested_event"
			return $condition
		}

		set genet_child		[dict get $event_info genet_child]

		# Duplications & speciations can be processed down
		set child_events [lsearch -all -inline -glob $events_list \{$genet_child,*\}*]

		if {[llength $child_events] == 1} {
			set should_be_dup	[ParseEvent [lindex $child_events 0]]
			if {[dict get $should_be_dup event_type] != "Dup"} {
				error "I think this must be a duplication. Code: 564"
				exit
			}
		} elseif {[llength $child_events] != 2} {
			error "I think there should be two events here: $child_events"
		}

		set events_to_test [concat $events_to_test $child_events]
		set events_to_test [lremove $events_to_test $tested_event]
	}

	return $condition
}


# proc IsTransferInOut {transfer_event parent_node child_node} {
# 	set final_subevent	[GetFinalSubevent $transfer_event]
# 	set parsed_transfer		[ParseTransferSubevent $final_subevent]
# 	set receiver_edge 		[dict get $parsed_transfer receiver_edge]
# 	# set receiver_nodes 		[split $receiver_edge \,]


# 	if {$receiver_edge == "$parent_node,$child_node"} {
# 		set transfer_type IN
# 	} else {
# 		set transfer_type OUT
# 	}

# 	return $transfer_type
# }




# proc PruneTransferEvents {events} {
# 	set pruned_events {}

# 	foreach event $events {
# 		set final_subevent [GetFinalSubevent $event]
# 		if {[IsSubeventTransfer $final_subevent] == FALSE} {
# 			lappend pruned_events $event
# 		}
# 	}
# 	return $pruned_events
# }

# proc transfer_subevent_parse {subevent} {
# 	global trans_type
# 	global gene_tree_xferred_child
# 	global donor_edge
# 	global receiver_edge

# 	set trans_type [string trim [string range $subevent [string first " " $subevent] end]]
# 	set subevent_split [split [string range $subevent 1 [string first " " $subevent]-2] \;]
# 	set gene_tree_xferred_child [lindex $subevent_split 0]
# 	set donor_edge [string range [lindex $subevent_split 1] 1 end-1]
# 	set receiver_edge [string range [lindex $subevent_split 2] 1 end-1]
# 	return
# }

# proc non_trans_subevent_parse {subevent} {
# 	global edge
# 	global child
# 	global parent
# 	global event_type

# 	set event_type [string trim [string range $subevent [string first " " $subevent] end]]
# 	set edge [string range $subevent 0 [string first " " $subevent]-1]
# 	set subevent_split [split [string range $edge 1 end-1] \,]
# 	set parent [lindex $subevent_split 0]
# 	set child [lindex $subevent_split 1]
# 	return
# }

# proc reduce_non_redundant_refine {redundant_events edge} {
# 	set non_redun_output {}

# 	## Put all the loss events at the end ##
# 	set temp_loss_list {}
# 	foreach event $redundant_events {
# 		set event_type [GetEventType $event]

# 		if {$event_type eq "Loss"} {
# 			set redundant_events [lremove $redundant_events $event]
# 			lappend temp_loss_list $event
# 		}
# 	}
# 	set redundant_events [concat $redundant_events $temp_loss_list]

# 	foreach event $redundant_events {
# 		set parsed_event	[ParseNontransSubevent $event]
# 		set edge 			[dict get $parsed_event edge]

# 		if {[lsearch -glob $non_redun_output *$edge*] == -1} {
# 			lappend non_redun_output $event
# 		}
# 	}
# 	return $non_redun_output
# }
