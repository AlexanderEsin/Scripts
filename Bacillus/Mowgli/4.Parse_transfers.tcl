#!/usr/local/bin/tclsh
#############
### Procs ###
#############

source ~/Documents/Scripts/General_utils.tcl
source ~/Documents/Scripts/Procs/mowgli_parsing.tcl

## Paths
set direct		/Users/aesin/Desktop/Bacillus/Mowgli/Mowgli_output
set out_dir		$direct/Raw_predictions
file mkdir		$out_dir/Logs


set transfer_cost		[lindex $argv 0]

set mowgli_output_dir	$direct/Output_$transfer_cost
set penalty_log_dir		$out_dir/Logs/$transfer_cost
file mkdir				$penalty_log_dir

set global_hgt_table	{}
set testedGroup			0
set OutAG_present		0
set global_OutAG		0
set scenario_2			0
set scenario_3			0

cd $mowgli_output_dir
set reconciliations [glob -type d *]

foreach reconciliation $reconciliations {

	# Log file for this gene family at this penalty
	set dir_log [open $penalty_log_dir/$reconciliation\.log w]

	# List of species tree nodes that are part of Anoxybacillus / Geobacillus branches
	set anoxy_geo_nodes [split [string trim [openfile $mowgli_output_dir/$reconciliation/coreNodes.tsv]] \n]
	
	multiputs $dir_log stdout "###########################################\nTree tested: $reconciliation"

	# Process the input mapping file to get a list of all the events
	set events_list [ParseMappingFile $mowgli_output_dir/$reconciliation/Fullmapping.mpr]
	if {[llength $events_list] == 0} {
		multiputs $dir_log stdout "Fullmapping.mpr file is empty. Mowgli run is incomplete / corrupted -- skipping $reconciliation directory"
		close $dir_log
		continue
	} else {
		incr testedGroup
	}

	# This will hold the output
	set HGT_table_out [list "Gene Family\tTransfer Index\tHGT Type\tDonor Edge\tReceiver Edge\tGeneT Parent\tGeneT Child\tHeritage\tTransfer Event"]
	
	# Holders
	set anoxy_geo_events	{}

	set transfer_AGOut	{}
	set transfer_AG2AG	{}
	set transfer_OutAG	{}

	set OutAG_done			{}
	set AG2AG_process		{}
	set AG2AG_done			{}

	set NonTrans_AG			{}
	set NonTrans_AG_process	{}
	set NonTrans_AG_done	{}

	set scenario			""
	set AG_root				FALSE
	set edge_heritage	[dict create "event" "heritage"]

	foreach node $anoxy_geo_nodes {

		set events [lsearch -all -inline -glob $events_list *,$node\)*]
		if {[llength $events] == 0} {
			continue
		}
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
					# Transfer from AG into non-AG branches
					lappend transfer_AGOut $event
				} elseif {[lsearch_any $anoxy_geo_nodes $donor_nodes] == TRUE} {
					# Transfer from AG into AG
					lappend transfer_AG2AG $event
				} else {
					# Transfer from non-AG branches into AG
					lappend transfer_OutAG $event
				}
			} else {
				lappend NonTrans_AG $event
			}
		}
	}

	set total_into_AG [expr [llength $transfer_AG2AG] + [llength $transfer_OutAG]]

	set explained_by_transfer 	0
	set transfer_counter		1

	foreach non_trans_AG_event $NonTrans_AG {
		set parsed_event		[ParseEvent $non_trans_AG_event]
		set event_parent_node	[dict get $parsed_event genet_parent]
		set edge				[dict get $parsed_event edge]

		# This edge is the tree root - it lies within AG. Will be reprocessed below
		if {[llength $event_parent_node] == 0} {
			lappend NonTrans_AG_process $non_trans_AG_event
			continue
		}

		set parent_child_events	[GetParentEvent $non_trans_AG_event $events_list]
		set parent_event		[dict get $parent_child_events parent_event]
		set parsed_parent		[ParseEvent $parent_event]
		
		# In v. rare situations the 'root' is actually a transfer event 'out'. As a result, the correct interpretation is that the child event (i.e. input) is the root. E.g. 7165 T=3	
		set root_flag			[IsEventRoot $parent_event]
		if {[dict get $parsed_parent type] eq "Trans" && $root_flag == TRUE && [IsTransIntoEdge $parent_event $edge] == FALSE} {
			puts $non_trans_AG_event
			# puts $dir_log "\t--> Event: $non_trans_AG_event has a parent that is both a Root and a Transfer Out:\n\t$parent_event"
			lappend NonTrans_AG_done "--> Event: $non_trans_AG_event has a parent that is both a Root and a Transfer Out:\n\t-->$parent_event"
			dict set edge_heritage $parent_event "Root"
			set AG_root TRUE
			incr scenario_3
			continue
		}

		# If the parent of this event (speciation, duplication, extant) is a transfer, it must be a transfer into AG
		if {[dict get $parsed_parent type] eq "Trans"} {
			
			if {[lsearch $transfer_OutAG $parent_event] != -1} {
				set transfer_event	[lindex [lsearch -all -inline $transfer_OutAG $parent_event] 0]
				set transfer_type	"OutAG"
			} elseif {[lsearch $transfer_AG2AG $parent_event] != -1} {
				set transfer_event	[lindex [lsearch -all -inline $transfer_AG2AG $parent_event] 0]
				set transfer_type	"AG2AG"
			} else {
				error "If this AG event has a Transfer parent event - it must be AG2AG or OutAG. Code 2222"
				exit
			}
			set parsed_transfer	[ParseEvent $transfer_event]
			set donor_edge		[dict get $parsed_transfer donor_edge]
			
			# Here we want to test whether the transfer is predicted into the most parsimonious branch. If not, find the most parsimonious transfer
			set true_receiver_event		[ReduceTransferBranch $non_trans_AG_event $transfer_AG2AG $events_list $dir_log $reconciliation]

			set parsed_receiver			[ParseEvent $true_receiver_event]
			set receiver_genet_parent	[dict get $parsed_receiver genet_parent]
			set receiver_genet_child	[dict get $parsed_receiver genet_child]
			set receiver_edge			[dict get $parsed_receiver edge]

			if {$transfer_type eq "OutAG"} {
				lappend OutAG_done	"$reconciliation\t$transfer_counter\tOutAG\t$donor_edge\t$receiver_edge\t$receiver_genet_parent\t$receiver_genet_child\tNA\t\|\|\t$parent_event"
				incr transfer_counter
			} else {
				lappend AG2AG_process	"HGT(AG2AG)\t$donor_edge\t$receiver_edge\t$receiver_genet_parent\t$receiver_genet_child\t\|\|\t$parent_event"
			}
			
			incr explained_by_transfer
			continue
		}

		set parent_child_node	[dict get $parsed_parent child]
		if {[lsearch $anoxy_geo_nodes $parent_child_node] == -1} {
			lappend NonTrans_AG_process $non_trans_AG_event
		}
	}

	if {$total_into_AG != $explained_by_transfer} {
		error "12345. Dir = $reconciliation"
		exit
	}

	# In certain cases, the predicted transfer will be into a deeper AG node, with
	# multiple losses futher on. A similar scenario (2) is considered in the context
	# of transfers above the common AG root below. Within AG such scenarious present
	# a non-parsimonious receptor node, and exist so that Mowgli can predict a
	# transfer out to another node from the "lost" taxa. These transfers from AG
	# "lost" taxa are unlikely to be informative, yet they will artificially
	# introduce more "ancient" HGT receipts than would otherwise be. So, for any
	# transfers into a non-terminal (not a single taxon) AG node need to be checked,
	# and refined down if necessary. This can be accomplished by checking whether
	# any of the two daughter edges post-transfer are lost entirely; if they are,
	# the new most-parsimonious receptor node becomes the non-lost edge, and the
	# process is repeated.

	# As we traverse the branches up, each edge can have a number of events
	# associated with it. Assuming the edges can be followed down to extant species
	# (which is always the case in our bottom up analysis), it is true that any
	# directly ancestral edge must then split into two other edges - UNLESS the edge
	# in question was horizontally transferred from another lineage. Due to the
	# ultrametric nature of the species tree - which is required for the model-based
	# assessment - but incorrect due to equivalent branch lengths assigned to the
	# tree (a quartet based method has no way to assign age-inferred branch
	# lengths). This arbitrary assignment of branch lengths is interpreted by mowgli
	# as the relative ancestral age, so in reality more broadly represented lineages
	# will typically have higher order branching and thus appear "more ancient" to
	# mowgli. When inferring direction of HGT, mowgli has to account for this
	# relative age discrepancy - and to make some of these gene trees work, the
	# program will actually predict a more ancestral transfer into a lineage,
	# followed by multiple losses (and perhaps transfers -out- to younger lineages).
	# This is not an issue, since the closest donor organism will be estimated by
	# parametric methods later; however, an HGT event into an ancestral node
	# followed by losses in all branches apart from the present Anoxy/Geobacillus
	# species should ALSO be considered as a transfer into the Anoxy/Geo taxa: eg
	# 1781 T=5

	# Now for each speciation that resulted in a topmost anoxy_geo node, we want to check whether the sister lineages have been lost


	
	# There can be multiple events at the root - this is because of duplications above root or transfers above root in the presence of a vertical homolog. The history of each event at the AG root is tested independetly

	foreach tested_event $NonTrans_AG_process {		
		set initial_event $tested_event
		set sister_fate		"Unknown"
		set scenario		""

		while {$sister_fate eq "Loss" || $sister_fate eq "Unknown"} {

			if {$sister_fate ne "Unknown"} {
				set tested_event	$parent_event
			}

			set parsed_event		[ParseEvent $tested_event]
			set event_parent_node	[dict get $parsed_event genet_parent]

			# If there is no parent node to this event, we are at the root
			if {[llength $event_parent_node] == 0} {
				set scenario 3
				break
			}

			set parent_child_events	[GetParentEvent $tested_event $events_list]
			set parent_event 		[dict get $parent_child_events parent_event]
			set parsed_parent		[ParseEvent $parent_event]

			# If parent type is a transfer, the origin of this branch is a transfer IN
			if {[dict get $parsed_parent type] eq "Trans"} {
				set scenario 2
				break
			}

			set penul_event		[dict get $parent_child_events penul_event]
			set sister_event	[GetSisterEvent $parent_event $penul_event $events_list]

			set sister_fate		[TestAllChildren $sister_event $transfer_AG2AG $events_list]
		}

		if {$scenario eq ""} {
			#puts $dir_log "\t--> Event: $initial_event has a Vertical history"
			lappend NonTrans_AG_done "--> Event: $initial_event has a Vertical history"
			dict set edge_heritage $initial_event "Vertical"
		# In scenario 2 the parent event is a TransIn
		} elseif {$scenario == 2} {
			set reduced_transfer		[ReduceTransferBranch $tested_event $AG2AG_transfer $events_list $dir_log $reconciliation]
			set parsed_reduced			[ParseEvent $reduced_transfer]

			set donor_edge				[dict get $parsed_parent donor_edge]	
			set receiver_edge			[dict get $parsed_reduced edge]
			set receiver_genet_parent	[dict get $parsed_reduced genet_parent]
			set receiver_genet_child	[dict get $parsed_reduced genet_child]
			
			lappend OutAG_done	"$reconciliation\t$transfer_counter\tOutAG\t$donor_edge\t$receiver_edge\t$receiver_genet_parent\t$receiver_genet_child\tNA\t\|\|\t$parent_event"
			incr transfer_counter
			incr total_into_AG

			#puts $dir_log "\t--> Event: $initial_event is an HGT above AG root with multiple losses (Scenario 2)"
			lappend NonTrans_AG_done "--> Event: $initial_event is an HGT above AG root with multiple losses (Scenario 2)"
			incr scenario_2
		} elseif {$scenario == 3} {
			#puts $dir_log "\t--> Event: $initial_event is part of the root (Scenario 3)"
			lappend NonTrans_AG_done "--> Event: $initial_event is part of the root (Scenario 3)"
			dict set edge_heritage $initial_event "Root"
			set AG_root TRUE
			incr scenario_3
		} else {
			error "There can be no other fates. Directory: $reconciliation. Code 9999"
			exit
		}
	}




	while {[llength $AG2AG_process] > 0} {
		set AG2AG_transfer [lindex $AG2AG_process 0]
		set delay_condition FALSE

		set info			[wsplit $AG2AG_transfer ||]
		set trans_event		[string trimleft [lindex $info end]]
		set parsed_event	[ParseEvent $trans_event]
		set test_edge		[dict get $parsed_event donor_edge]

		# This is needed to correctly format the output table
		set edge_info		[string trim [lindex $info 0]]
		set edge_output	[join [lrange [split $edge_info \t] 1 end] \t]

		# The transfer out can be the root, but in that case it's already been added to the heritage dictionary
		if {[IsEventRoot $trans_event] == TRUE} {
			set AG2AG_heritage	[dict get $edge_heritage $trans_event]

			lappend AG2AG_done "$reconciliation\t$transfer_counter\tAG2AG\t$edge_output\t$AG2AG_heritage\t\|\|\t$trans_event"
			set AG2AG_process [lremove $AG2AG_process $AG2AG_transfer]
			incr transfer_counter
			continue
		}

		# Get the parental non-TransOut event to this AG2AG transfer
		set parent_not_TransOut	[SkipAllTransOut $trans_event $test_edge $events_list "up"]
		set parent_event		[dict get $parent_not_TransOut parent_event]

		# If the immediate parent is not TransIn then then it can only be a Spec/Dup
		while 1 {		
		
			set heritage_known	[dict exists $edge_heritage $parent_event]
			if {$heritage_known == 1} {
				set AG2AG_heritage	[dict get $edge_heritage $parent_event]
				lappend AG2AG_done "$reconciliation\t$transfer_counter\tAG2AG\t$edge_output\t$AG2AG_heritage\t\|\|\t$trans_event"
				incr transfer_counter
				break
			}

			set parsed_parent	[ParseEvent $parent_event]
			# If the event type is a transfer - it must be a TransIn (i.e. an OutAG)
			if {[dict get $parsed_parent type] eq "Trans"} {

				# Find the corresponding OutAG entry
				set parsed_test		[ParseEvent $parent_event]
				set genet_parent	[dict get $parsed_test genet_parent]
				set genet_child		[dict get $parsed_test genet_child]

				## If it's a TransIn - it could be from out (we can find it in $OutAG_done) or it can be from AG also (look in AG2AG done). It's possible this AG2AG is being tested before the one its descended from, in which case we simply move this AG2AG to the back of the processing list and continue

				# Try to find it in OutAG
				set match_OutAG [lsearch -all -inline -glob $OutAG_done "*||\t$parent_event*"]
				if {[llength $match_OutAG] == 1} {
					set match		[lindex $match_OutAG 0]
					set AG2AG_heritage	[lindex [split $match \t] 1]
				# Try to find it in AG2AG
				} else {
					set match_AG2AG [lsearch -all -inline -glob $AG2AG_done "*||\t$parent_event*"]
					if {[llength $match_AG2AG] == 1} {
						set match		[lindex $match_AG2AG 0]
						set AG2AG_heritage	[lindex [lindex [wsplit $match ||] 0] 7]
					} else {
						set delay_condition TRUE
						break
					}
				}

				lappend AG2AG_done "$reconciliation\t$transfer_counter\tAG2AG\t$edge_output\t$AG2AG_heritage\t\|\|\t$trans_event"
				incr transfer_counter
				break
			}

			set parent_child_node	[dict get $parsed_parent child]
			set parent_event	[dict get [GetParentEvent $parent_event $events_list] parent_event]
		}

		set AG2AG_process [lremove $AG2AG_process $AG2AG_transfer]
		if {$delay_condition == TRUE} {
			lappend AG2AG_process $AG2AG_transfer
		}
	}

	set total_finished [expr [llength $AG2AG_done] + [llength $OutAG_done]]
	if {$total_into_AG != $total_finished} {
		error "P1GS"
		exit
	}

	# Process the HGT output lines for easy table read
	set all_HGT_done	[concat $OutAG_done $AG2AG_done]
	set formatted_out	[FormatHgtOutput $all_HGT_done]
	set HGT_table_out	[concat $HGT_table_out $formatted_out]


	multiputs stdout $dir_log [join [concat $HGT_table_out $NonTrans_AG_done] \n]
	set global_OutAG	[expr $global_OutAG + [llength $OutAG_done]]

	##
	## Per gene family output
	##

	if {[llength $OutAG_done] != 0} {
		set HGT_present 1
		incr OutAG_present
		lappend global_hgt_table "$reconciliation\t1\t$AG_root"
	} else {
		set HGT_present 0
		lappend global_hgt_table "$reconciliation\t0\t$AG_root"
	}

	# Write out the file containing all the donor and receiver edges for the transfers INTO geobacillus from an external source ##
	set out [open $mowgli_output_dir/$reconciliation/Parsed_events.tsv w]
	puts $out [join $HGT_table_out \n]
	close $out

	# Write out the log file
	close $dir_log
}
		
set global_hgt_table [concat [list "Gene Family\tHGT Present\tIs AG root"] [lsort -dictionary $global_hgt_table]]
set out_table_filename $out_dir/T$transfer_cost\_D2_L1_results.txt
set out_stats_filename $out_dir/T$transfer_cost\_D2_L1_stats.txt

set out [open $out_table_filename w]
puts $out [join $global_hgt_table \n]
close $out

set stats [open $out_stats_filename w]
multiputs stdout $stats		"Total number of gene families tested: $testedGroup\nFamilies containing at least 1 OutAG: $OutAG_present\n\nTotal number of OutAG events: $global_OutAG\nTotal number of times where transfer into AG predicted above root (S2): $scenario_2\nTotal number of times where root is in AG (S3): $scenario_3\n"
close $stats	



