## Procs for Mowgli output parsing ##

## //// THIS NEEDS TO BE TIDIED 21/7/2017 - some obsolete functions need to be deleted //// ##

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

proc IsSubeventTransfer {subevent} {
	set is_transfer FALSE
	if {[regexp {;} $subevent] == 1} {
		set is_transfer TRUE
	}
	return $is_transfer
}

proc IsEventRoot {event} {
	set parsed_event	[ParseEvent $event]
	set genet_parent	[dict get $parsed_event genet_parent]
	if {[llength $genet_parent] == 0} {
		return TRUE
	} else {
		return FALSE
	}
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

# Skip all transfer events out of the lineage being tested in a linear way. If going up (high order / more ancient branches) then we can stop at a transfer event coming IN. If we are moving down, there cannot be another transfer IN event along the same path.
proc SkipAllTransOut {trans_event test_edge events_list {direction down}} {
	set next_event $trans_event

	set parsed_event [ParseEvent $next_event]
	set type [dict get $parsed_event type]

	# This loop has to run at least once because the event provided must be a transfer OUT #
	set root_flag FALSE
	while {$type eq "Trans" && [IsTransIntoEdge $next_event $test_edge] == FALSE} {
		set previous_event $next_event

		if {$root_flag == TRUE} {
			break
		}

		# For every transfer up - there can only be one parent
		if {$direction eq "up"} {
			set genet_parent	[dict get $parsed_event genet_parent]
			set next_event		[lsearch -all -inline -glob $events_list "\{*,$genet_parent\}*"]
			# We might skip up to the root
			if {[llength $next_event] == 0} {
				set next_event		[lsearch -all -inline -glob $events_list "\(*,$genet_parent\)*"]
				set root_flag		TRUE
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
		set root_flag		TRUE
		set parent_event	[lsearch -all -inline -glob $events_list "\( ,$child_genet_parent\)\t*"]
	} else {
		set root_flag		FALSE
	}

	# If we still can't find a parent event - error.
	if {[llength $parent_event] == 0} {
		error "Expecting one parent event. Code: D0GS"
	}

	set parent_event		[lindex $parent_event 0]
	set parsed_parent		[ParseEvent $parent_event]

	if {[dict get $parsed_parent type] eq "Trans" && [IsTransIntoEdge $parent_event $child_edge] == FALSE} {
		# In v. rare situations the 'root' is actually a transfer event 'out'. E.g. 2085 T=5
		if {$root_flag == TRUE} {
			set parent_child_events [dict create parent_event $parent_event]
			dict set parent_child_events penul_event $child_event
			dict set parent_child_events child_event $child_event
		} else {
			set parent_child_events [SkipAllTransOut $parent_event $child_edge $events_list "up"]
			dict set parent_child_events child_event $child_event
		}
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

proc TestAllChildren {start_event ag2agTrans_list events_list} {
	set events_to_test [list $start_event]

	set condition "Loss"
	while {[llength $events_to_test] != 0} {

		set tested_event	[lindex $events_to_test 0]
		# puts "---> Testing event: $tested_event"

		set event_info		[ParseEvent $tested_event]
		set type			[dict get $event_info type]

		if {$type eq "Trans"} {

			# If while reducing the receiver branch we find an AG2AG transfer. We don't want to reduce past this point
			if {[lsearch $ag2agTrans_list $tested_event] != -1} {
				set condition "AG2AG Trans: $tested_event"
				return $condition
			}

			set test_edge			[dict get $event_info donor_edge]
			set child_parent_events	[SkipAllTransOut $tested_event $test_edge $events_list "down"]
			set next_event			[dict get $child_parent_events child_event]
			set event_info			[ParseEvent $next_event]
			# puts "SKIPPED TRANSFER DOWN: $next_event"
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

proc ReduceTransferBranch {top_event ag2agTrans_list events_list dir_log out_dir} {
	while 1 {
		set parsed_event	[ParseEvent $top_event]
		set genet_child		[dict get $parsed_event genet_child]

		if {[dict get $parsed_event event_type] eq "Extant"} {
			return $top_event
			break
		}

		# Find the children events of the branch that received the transfer
		set child_events	[lsearch -all -inline -glob $events_list "\{$genet_child\,*"]

		foreach child_event $child_events {
			set fate	[TestAllChildren $child_event $ag2agTrans_list $events_list]
			# puts $fate
			if {$fate eq "Loss"} {
				set child_events	[lremove $child_events $child_event]
			}
		}

		# If there are still two child events - that means both branches down have extant taxa, if there is only one - there has been a total loss in one branch.
		if {[llength $child_events] == 2} {
			return $top_event
			break
		} elseif {[llength $child_events] == 1} {
			set next_event [lindex $child_events 0]

			# If a sister has been lost, it's possible that the next event in the remaining sister is a TransOut - we need to reduce this to find the next Spec/Dup/Ext (Loss?) event
			if {[dict get [ParseEvent $next_event] type] eq "Trans"} {
				set parsed_transfer	[ParseEvent $next_event]
				set donor_edge		[dict get $parsed_transfer donor_edge]
				set next_event		[dict get [SkipAllTransOut $next_event $donor_edge $events_list "down"] child_event]
			}

			# multiputs stdout $dir_log "One child of TransIn-derived event:\n\t$top_event\nhas been lost. Continuing down to: $next_event"
			set top_event $next_event
		} else {
			error "There should be some children left. Code 5555"
			exit
		}
	}	
}

proc FormatHgtOutput {HGT_events} {
	set formatted_list {}
	foreach HGT_event $HGT_events {
		set split_line		[wsplit $HGT_event ||]
		# Space separate the transfer event
		set format_trans	[join [split [string trim [lindex $split_line 1]] \t] " "]
		set other_cols		[split [string trim [lindex $split_line 0]] \t]
		lappend other_cols	$format_trans
		set HGT_formatted	[join $other_cols \t]
		lappend formatted_list $HGT_formatted
	}
	return $formatted_list
}
