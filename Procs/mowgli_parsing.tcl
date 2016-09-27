## Procs for Mowgli output parsing ##

proc parse_mapping_file {mapping_file} {
	global events_list
	global empty_map_switch

	set empty_map_switch 0
	## Read in the Fullmapping file ##
	set recon_map [string trim [openfile $mapping_file]]
	if {[string length $recon_map] == 0} {
		puts stderr "ERROR: mapping file is empty"
		set empty_map_switch 1
		return
	}
	## Parse the input data to remove the header lines, any extra tabs and newlines ##
	regsub -all {\n\n} $recon_map "\n" recon_map
	regsub -all {\t\t} $recon_map "\t" recon_map
	## Events_list is a list of all the individual evolutionary events as given by mowgli ##
	set events_list [lrange [split $recon_map \n] 2 end]
	return $events_list
}

proc get_transfer_events {events} {
	global transfer_events

	set transfer_events {}
	foreach event $events {
		if {[lsearch -glob $event *Tran*] > 0} {
			lappend transfer_events $event
		}
	}
	return $transfer_events
}

proc get_final_subevent {event} {
	global final_subevent

	set event_data [split [string trim $event] \t]
	set final_subevent [lindex $event_data end]
	return $final_subevent
}

proc transfer_subevent_parse {subevent} {
	global trans_type
	global gene_tree_xferred_child
	global donor_edge
	global receiver_edge

	set trans_type [string trim [string range $subevent [string first " " $subevent] end]]
	set subevent_split [split [string range $subevent 1 [string first " " $subevent]-2] \;]
	set gene_tree_xferred_child [lindex $subevent_split 0]
	set donor_edge [string range [lindex $subevent_split 1] 1 end-1]
	set receiver_edge [string range [lindex $subevent_split 2] 1 end-1]
	return
}

proc non_trans_subevent_parse {subevent} {
	global edge
	global child
	global parent
	global event_type

	set event_type [string trim [string range $subevent [string first " " $subevent] end]]
	set edge [string range $subevent 0 [string first " " $subevent]-1]
	set subevent_split [split [string range $edge 1 end-1] \,]
	set parent [lindex $subevent_split 0]
	set child [lindex $subevent_split 1]
	return
}

proc is_subevent_transfer {subevent} {
	if {[regexp {;} $subevent] == 1} {
		set transfer TRUE
	} else {
		set transfer FALSE
	}
	return $transfer
}

proc reduce_non_redundant_edges {redundant_events} {
	set non_redun_output {}
	global edge

	foreach event $redundant_events {
		set final_subevent [get_final_subevent $event]
		non_trans_subevent_parse $final_subevent

		if {[lsearch -glob $non_redun_output *$edge*] == -1} {
			lappend non_redun_output $event
		}
	}
	return $non_redun_output
}

proc reduce_non_redundant_refine {redundant_events} {
	set non_redun_output {}
	global edge

	## Put all the loss events at the end ##
	set temp_loss_list {}
	foreach event $redundant_events {
		set event_type [string trim [string range $event [string first " " $event] end]]

		if {$event_type eq "Loss"} {
			set redundant_events [lremove $redundant_events $event]
			lappend temp_loss_list $event
		}
	}
	set redundant_events [concat $redundant_events $temp_loss_list]

	foreach event $redundant_events {
		non_trans_subevent_parse $event

		if {[lsearch -glob $non_redun_output *$edge*] == -1} {
			lappend non_redun_output $event
		}
	}
	return $non_redun_output
}

proc reduce_non_redundant_children {redundant_events} {
	set non_redun_output {}
	global edge

	foreach event $redundant_events {
		set final_subevent [get_final_subevent $event]
		non_trans_subevent_parse $final_subevent

		if {[lsearch -glob $non_redun_output *$edge\ *] == -1} {
			lappend non_redun_output $event
		}
	}
	return $non_redun_output
}