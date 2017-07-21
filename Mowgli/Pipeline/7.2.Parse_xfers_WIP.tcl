#!/usr/local/bin/tclsh
#############
### Procs ###
#############

source ~/Documents/Scripts/General_utils.tcl
source ~/Documents/Scripts/Procs/mowgli_parsing.tcl

set transfer_cost 5

## Paths
set direct /Users/aesin/Desktop/Mowgli/Mowgli_outputs

## Make the two output_directories for the results/stats files
set output_path	$direct/../Raw_predictions
file mkdir $output_path/Logs

## Regex patterns
set pata {\(+?[0-9]+?\,+?}
set patb {\)+?}


set mowgli_output_dir	$direct/Mow_test_out_t$transfer_cost
set penalty_log_dir		$output_path/Logs/$transfer_cost
file mkdir				$penalty_log_dir

set global_hgt_table {}

####################
### Set counters ###
####################

set all_tested 0
set vertical 0
set horizontal 0
set scenario_1_hgt 0
set scenario_2_hgt 0
set scenario_3_hgt 0
# Scenario 4 is a combination of 2 & 3 - all sisters are lost UP TO the root #
set scenario_4_hgt 0
set ag2ag_transfer_counter 0

##############################################################################

cd $mowgli_output_dir
set output_directories [glob -type d *]
set out_dir 1289
set output_directories [lrange $output_directories 0 99]


//// THIS ALL SEEMS TO WORK. NEED TO TEST FURTHER AND COMMENT. 21/7/2017 ////





# foreach out_dir $output_directories {

	####################################
	### Empty variables to be filled ###
	####################################
	set HGT_table_out "$out_dir"

	set scenario ""

	## Log file for this gene family at this penalty
	set dir_log [open $penalty_log_dir/$out_dir\.log w]

	cd $mowgli_output_dir/$out_dir
	multiputs $dir_log stdout "\n###########################################\nTree tested: $out_dir"

	## RUN the NODE identifier script here to get a list of all the nodes that correspond to anoxy/geobacillus in the species tree ##
	## Output file MUST be called Anoxy_geo_nodes.tsv ##
	set anoxy_geo_nodes [split [string trim [openfile Anoxy_geo_nodes.tsv]] \n]

	# Process the input mapping file to get a list of all the events
	set events_list [ParseMappingFile Fullmapping.mpr]
	if {$events_list == 0} {
		multiputs $dir_log stdout "Fullmapping.mpr file is empty. Mowgli run is incomplete / corrupted -- skipping $out_dir directory"
		exit
	}
	
	# Holders
	set anoxy_geo_events	{}
	set AG_present_nodes	{}

	set transfer_out_of_AG	{}
	set transfer_AG2AG		{}
	set transfer_into_AG	{}
	set non_trans_AG_events	{}
	set events_from_out		{}

	# 
	foreach node $anoxy_geo_nodes {

		set events [lsearch -all -inline -glob $events_list *,$node\)*]
		if {[llength $events] == 0} {
			continue
		}
		lappend AG_present_nodes $node

		foreach event $events {
			# In transfer events between AG nodes, the event will be picked up twice - once as the donor node, once as a receiver node. Exclude double entries.
			if {[lsearch $anoxy_geo_events $event] > -1 } {
				continue
			}

			lappend anoxy_geo_events $event
			set parsed_event	[ParseEvent $event]

			# If it's a transfer, make sure that the the receiver node is an AG node - otherwise can pick up AG -> other transfers.
			if {[dict get $parsed_event type] eq "Trans"} {
				set receiver_edge	[string range [dict get $parsed_event receiver_edge] 1 end-1]
				set donor_edge		[string range [dict get $parsed_event donor_edge] 1 end-1]

				set receiver_child	[lindex [split $receiver_edge \,] 1]
				set donor_nodes		[split $donor_edge \,]

				if {[lsearch $anoxy_geo_nodes $receiver_child] == -1} {
					lappend transfer_out_of_AG $event
				} elseif {[lsearch_any $anoxy_geo_nodes $donor_nodes] == TRUE} {
					lappend transfer_AG2AG $event
				} else {
					lappend transfer_into_AG $event
				}
			} else {
				lappend non_trans_AG_events $event
			}
		}
	}

	set total_into_AG [expr [llength $transfer_AG2AG] + [llength $transfer_into_AG]]

	puts [llength $transfer_out_of_AG]
	puts [llength $transfer_AG2AG]
	puts [llength $transfer_into_AG]
	puts [llength $non_trans_AG_events]

	set explained_by_transfer 0
	foreach non_trans_AG_event $non_trans_AG_events {
		set parsed_event	[ParseEvent $non_trans_AG_event]
		set event_parent_node	[dict get $parsed_event genet_parent]

		# This edge is the root - it lies within AG. Will be reprocessed below
		if {[llength $event_parent_node] == 0} {
			lappend events_from_out $non_trans_AG_event
			continue
		}

		set parent_child_events	[GetParentEvent $non_trans_AG_event $events_list]
		set parent_event		[dict get $parent_child_events parent_event]

		set parsed_parent		[ParseEvent $parent_event]
		if {[dict get $parsed_parent type] eq "Trans"} {
			set non_trans_AG_events [lremove $non_trans_AG_events non_trans_AG_event]
			incr explained_by_transfer
			continue
		}

		set parent_child_node	[dict get $parsed_parent child]
		if {[lsearch $anoxy_geo_nodes $parent_child_node] == -1} {
			lappend events_from_out $non_trans_AG_event
		}
	}

	if {$total_into_AG != $explained_by_transfer} {
		error "12345. Dir = $out_dir"
		exit
	}

	puts $explained_by_transfer
	puts $events_from_out


	# 		## In certain cases, the predicted transfer will be into a deeper AG node, with multiple losses futher on. A similar scenario (2) is considered in the context of transfers above the common AG root below. Within AG such scenarious present a non-parsimonious receptor node, and exist so that Mowgli can predict a transfer out to another node from the "lost" taxa. These transfers from AG "lost" taxa are unlikely to be informative, yet they will artificially introduce more "ancient" HGT receipts than would otherwise be. So, for any transfers into a non-terminal (not a single taxon) AG node need to be checked, and refined down if necessary. This can be accomplished by checking whether any of the two daughter edges post-transfer are lost entirely; if they are, the new most-parsimonious receptor node becomes the non-lost edge, and the process is repeated. ##

	## As we traverse the branches up, each edge can have a number of events associated with it.
	## Assuming the edges can be followed down to extant species (which is always the case in our bottom up analysis), it is true that any directly ancestral edge must then split into two other edges - UNLESS the edge in question was horizontally transferred from another lineage.
	## Due to the ultrametric nature of the species tree - which is required for the model-based assessment - but incorrect due to equivalent branch lengths assigned to the tree (a quartet based method has no way to assign age-inferred branch lengths). 
	## This arbitrary assignment of branch lengths is interpreted by mowgli as the relative ancestral age, so in reality more broadly represented lineages will typically have higher order branching and thus appear "more ancient" to mowgli. 
	## When inferring direction of HGT, mowgli has to account for this relative age discrepancy - and to make some of these gene trees work, the program will actually predict a more ancestral transfer into a lineage, followed by multiple losses (and perhaps transfers -out- to younger lineages). 
	## This is not an issue, since the closest donor organism will be estimated by parametric methods later; however, an HGT event into an ancestral node followed by losses in all branches apart from the present Anoxy/Geobacillus species should ALSO be considered as a transfer into the Anoxy/Geo taxa ##
	## Eg 1781 T=5

	# Now for each speciation that resulted in a topmost anoxy_geo node, we want to check whether the sister lineages have been lost ##

	if {[llength $events_from_out] == 1} {

		set tested_event	[lindex $events_from_out 0]
		set sister_fate		"Unknown"

		while {$sister_fate eq "Loss" || $sister_fate eq "Unknown"} {

			if {$sister_fate ne "Unknown"} {
				set tested_event	$parent_event
			}

			set parsed_event		[ParseEvent $tested_event]
			set event_parent_node	[dict get $parsed_event genet_parent]

			# If there is no parent node to this event, we are at the root
			if {[llength $event_parent_node] == 0} {
				puts "Early Break == IsRoot: 3"
				set scenario 3
				break
			}

			set parent_child_events	[GetParentEvent $tested_event $events_list]
			set parent_event 		[dict get $parent_child_events parent_event]
			set parsed_parent		[ParseEvent $parent_event]

			# If parent type is a transfer, the origin of this branch is a transfer IN
			if {[dict get $parsed_parent type] eq "Trans"} {
				puts "Early break == Transfer above edge with multiloss: 2"
				set scenario 2
				break
			}

			set penul_event		[dict get $parent_child_events penul_event]
			set sister_event	[GetSisterEvent $parent_event $penul_event $events_list]

			set sister_fate		[TestAllChildren $sister_event $events_list]
			puts "Current fate: $sister_fate"
		}
	}





	if {[llength $HGT_table_out] == 1} {
		lappend HGT_table_out "No transfers into Anoxy/Geobacillus"
		lappend global_hgt_table "$out_dir\t0\t$scenario"
		incr vertical
	} else {
		# lappend HGT_table_out "\n[join $ag2ag_event_l \n]"
		lappend global_hgt_table "$out_dir\t1\t$scenario"
		incr horizontal
	}
	# Write out the file containing all the donor and receiver edges for the transfers INTO geobacillus from an external source ##
	set out [open Transfers_in.tsv w]
	puts $out [join $HGT_table_out \n]
	close $out

	multiputs $dir_log stdout [join $HGT_table_out \n]

	## Write out the log file 
	close $dir_log

	incr all_tested
}
		

	# set global_hgt_table [lsort -dictionary $global_hgt_table]


	# set out_table_filename $output_path/T$transfer_cost\_D2_L1_results.txt
	# set out_stats_filename $output_path/T$transfer_cost\_D2_L1_stats.txt


	# set global_hgt_table [join $global_hgt_table \n]
	# set out [open $out_table_filename w]
	# puts $out $global_hgt_table
	# close $out

	# set stats [open $out_stats_filename w]
	# puts $stats "Total number of gene families tested: $all_tested\nAG vertical: $vertical\nAG horizontal: $horizontal\n\nScenario 1: $scenario_1_hgt\nScenario 2: $scenario_2_hgt\nScenario 3: $scenario_3_hgt\nScenario 4: $scenario_4_hgt"
	# close $stats



