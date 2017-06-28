#!/usr/local/bin/tclsh

###########################################################################
## Procs ##

source ~/Dropbox/Scripts/General_utils.tcl
source ~/Dropbox/Scripts/Procs/mowgli_parsing.tcl

###########################################################################

set transfer_cost 5
set ignore_inside_only_switch 0

####################################

## Regex patterns ##
set pata {\(+?[0-9]+?\,+?}
set patb {\)+?}

####################################

## Set script-wide variables ##
set direct /Users/aesin/Desktop/Mowgli/Mowgli_outputs
set mowgli_output_dir $direct/Mow_test_out_t$transfer_cost
set global_hgt_table {}

#######################################

## Make the two output_directories for the results/stats files ##
file mkdir $direct/Results_Full
file mkdir $direct/Results_NO_IG

#######################################

## Set counters for the different events ##
set all_tested 0
set vertical 0
set horizontal 0
set scenario_1_hgt 0
set scenario_2_hgt 0
set scenario_3_hgt 0
# Scenario 4 is a combination of 2 & 3 - all sisters are lost UP TO the root #
set scenario_4_hgt 0

#######################################

## Read in the file containing a list of IG-only trees ##
if {$ignore_inside_only_switch == 1} {
	cd $direct/..
	set IG_trees_list [split [string trim [openfile Trees_in_group_only.tsv]] \n]
}


#######################################

cd $mowgli_output_dir
set output_directories [glob -type d *]

foreach out_dir $output_directories {
	if {$ignore_inside_only_switch == 1} {
		if {[lsearch $IG_trees_list $out_dir] != -1} {
			puts "This tree: $out_dir contains only IG taxa - skipping"
			continue
		}
	}

	set HGT_table_out "$out_dir"

	cd $mowgli_output_dir/$out_dir
	puts "\n###########################################\nTree tested: $out_dir"

	## RUN the NODE identifier script here to get a list of all the nodes that correspond to anoxy/geobacillus in the species tree ##
	## Output file MUST be called Anoxy_geo_nodes.tsv ##
	set anoxy_geo_nodes [split [string trim [openfile Anoxy_geo_nodes.tsv]] \n]

	## Process the input mapping file to get a list of all the events #
	parse_mapping_file Fullmapping.mpr
	if {$empty_map_switch == 1} {
		puts "Run incomplete -- skipping $out_dir directory"
		continue
	}
	puts "Total number of events: [llength $events_list]"

	## Isolate just the transfer_events ##
	get_transfer_events $events_list
	puts "Number of transfer events: [llength $transfer_events]\n"

	###########################################

	## First we need to map the speciation history of any Anoxy/Geo taxa ##
	set speciation_from_outside {}
	set dups_form_outside {}
	set transfer_from_outside {}
	set transfer_from_anoxygeo {}
	set anoxygeo_present_nodes {}
	set root_multi_loss 0
	set scenario ""
	set fate ""

	## Each node that is present in the tree will have a map entry, thus any anoxy_geo nodes (or tips) that don't appear in the events are not present in the tree ##
	foreach node $anoxy_geo_nodes {
		set events [lsearch -all -inline -glob $events_list *,$node\)*]
		if {[llength $events] != 0} {
			lappend anoxygeo_present_nodes $node
		}
	}

	#puts $anoxygeo_present_nodes
	foreach node $anoxygeo_present_nodes {
		## Find all events for each anoxy_geo node ##
		set events [lsearch -all -inline -glob $events_list *,$node\)*]
		#puts $events
	
		if {[llength $events] != 0} {
			foreach event $events {

				## By looking for just the final_subevent we are excluding those cases where there has been a transfer out of a lineage, followed by a loss in that lineage ##
				## However, since a transfer out of a branch of Geobacillus that has then been lost is not interesting (right now), as we are looking only for transfers in, we can ignore that scenario for the time being ##

				set final_subevent [get_final_subevent $event]

				## Only transfers contain a semicolon ##
				if {[is_subevent_transfer $final_subevent] == TRUE} {

					transfer_subevent_parse $final_subevent
					set donor_edge_ends [split $donor_edge \,]
					## For the transfer, is this node part of the donor? If so, the transfer is OUT of Geobacillus ##
					if {[lsearch $anoxygeo_present_nodes [lindex $donor_edge_ends 0]] != -1 || [lsearch $anoxygeo_present_nodes [lindex $donor_edge_ends 1]] != -1} {
						## This node has had a transfer from within the anoxy_geo clade ##
						lappend transfer_from_anoxygeo $event
					} else {
						## This node has had a transfer from outside the anoxy_geo clade ##
						lappend transfer_from_outside $event
					}

				## Everything else is either a speciation, is extant, is a duplication, or a loss ##
				} else {
					non_trans_subevent_parse $final_subevent

					if {$event_type eq "Loss"} {
						continue
					} elseif {$event_type eq "NoEvent"} {
						puts "Not sure that this should happen - code: A1"
						exit 2
					} elseif {$event_type eq "Dup"} {
						if {[lsearch $anoxygeo_present_nodes $parent] != -1} {
							## The parent of this duplication is an anoxy_geo_node ##
							continue
						} else {
							## The parent of this duplication is not an anoxy_geo_node ##
							lappend dups_form_outside $event
						}
					} else {
						if {[lsearch $anoxygeo_present_nodes $parent] != -1} {
							## The parent of this speciation or extant taxon is an anoxy_geo_node ##
							continue
						} else {
							## The parent of this speciation or extant taxon is not an anoxy_geo_node ##
							lappend speciation_from_outside $event
						}
					}
				}
			}
		}
	}



	## Collapse redundant duplication edges - i.e. an edge that had multiple duplication events ##
	set non_redundant_dups_from_outside [reduce_non_redundant_edges $dups_form_outside]
	puts "Number of duplications of an edge not descendent from a present node: [llength $dups_form_outside]"
	puts "Reduced to non-redundant: [llength $non_redundant_dups_from_outside]"

	## Collapse redundant speciation edges - i.e. an edge that had multiple duplication, then as a results, had multiple speciation events ##
	set non_redundant_spec_from_outside [reduce_non_redundant_edges $speciation_from_outside]
	puts "Number of speciations of an edge not descendent from a present node: [llength $speciation_from_outside]"
	puts "Reduced to non-redundant: [llength $non_redundant_spec_from_outside]"

	## Combine speciation events from an external source with duplication events from an external source ##
	set events_from_outside [concat $non_redundant_spec_from_outside $non_redundant_dups_from_outside]
	puts "Number of duplication/speciation events from external sources: [llength $events_from_outside]"

	## Collapse all the connected speciation and duplication events ##
	puts "\nAttempting to collapse any redundant speciation and duplication events: [llength $events_from_outside]"
	if {[llength $events_from_outside] > 1} {
		foreach external_speciation $events_from_outside {
			set final_subevent [get_final_subevent $external_speciation]
			non_trans_subevent_parse $final_subevent
			set external_speciation_edge $edge

			if {$event_type eq "Dup"} {
				if {[lsearch -glob $events_from_outside "*$edge\ Spec*"] != -1 || [lsearch -glob $events_from_outside "*$edge\ Ext*"] != -1} {
					puts "--> Edge $external_speciation_edge duplicated then speciated / remained extant"
					set events_from_outside [lremove $events_from_outside $external_speciation]
				}
			}
		}
	}


	## let's take a shortcut here ##
	## If we find a speciation event that the product of which is the same edge as from an external transfer - we can summarise those to just the transfer ##
	puts "\nAttempting to assign transfers to the speciation/duplication events ..."
	foreach external_event $events_from_outside {
		set final_subevent [get_final_subevent $external_event]
		non_trans_subevent_parse $final_subevent
		set external_speciation_edge $edge

		if {[lsearch -glob $transfer_from_outside *\;$edge*] > -1} {			
			set transfer_hit [lsearch -inline -glob $transfer_from_outside *\;$edge*]
			set transfer_final_subevent [get_final_subevent $transfer_hit]
			transfer_subevent_parse $transfer_final_subevent

			##########################

			## In certain cases, the predicted transfer will be into a deeper AG node, with multiple losses futher on. A similar scenario (2) is considered in the context of transfers above the common AG root beow. Within AG such scenarious present a non-parsimonious receptor node, and exist so that Mowgli can predict a transfer out to another node from the "lost" taxa. These transfers from AG "lost" taxa are unlikely to be informative, yet they will artificially introduce more "ancient" HGT receipts than would otherwise be. So, for any transfers into a non-terminal (not a single taxon) AG node need to be checked, and refined down if necessary. This can be accomplished by checking whether any of the two daughter edges post-transfer are lost entirely; if they are, the new most-parsimonious receptor node becomes the non-lost edge, and the process is repeated. ##

			set test_the_edge $receiver_edge

			while 1 {

				## Obviously if the transfer is into a single taxon - there can be no lower node ##
				if {$event_type eq "Extant"} {
					break
				}

				set transfered_edge_child_node [lindex [split $test_the_edge \,] 1]
				set children [lsearch -all -inline -glob $events_list *\($transfered_edge_child_node,*]
				set branch_fates {}

				foreach child_event $children {
					set final_subevent [get_final_subevent $child_event]
					if {[is_subevent_transfer $final_subevent] == TRUE} {
						continue
					} else {
						lappend branch_fates $final_subevent
					}
				}

				## If there are duplications, for each duplication picked up - there will be two further speciations - we can reduce this down to a single entry for the fate ##
				set branch_fates [reduce_non_redundant_refine $branch_fates]
				## Sanity check ##
				if {[llength $branch_fates] != 2} {
					puts "Error: not two fate events. Tree: $out_dir\t$receiver_edge A1"; exit
				}

				## Remove any branches that have a loss ##
				foreach fate $branch_fates {
					non_trans_subevent_parse $fate
					if {$event_type eq "Loss"} {
						set branch_fates [lremove $branch_fates $fate]
					}
				}

				## Check whether both speciated, or one was lost ##
				if {[llength $branch_fates] == 2} {
					puts "This is the final edge: $test_the_edge"
					set receiver_edge $test_the_edge
					break
				} elseif {[llength $branch_fates] == 0} {
					puts "Error: two fates resulted in two losses. Tree: $out_dir A2"; exit
				} else {
					non_trans_subevent_parse [lindex $branch_fates 0]
					## If of two fates, one is loss and the other is extant, there can be no further down-sampling of the tree. We take the transfer to come into the extant taxon ##
					if {$event_type eq "Extant"} {
						set receiver_edge [string range $edge 1 end-1]
						puts "This is the final edge: $receiver_edge"
						break
					} else {
						puts "One branch was lost, continuing down..."
						set test_the_edge [string range $edge 1 end-1]
					}		
				}
			}

			##########################

			puts "This edge $external_speciation_edge was transferred horizontally from an external source"
			lappend HGT_table_out "HGT: $donor_edge\t$receiver_edge"
			incr scenario_1_hgt
			set scenario 1

			set events_from_outside [lremove $events_from_outside $external_event]

			set xfer_to_rm [lsearch -inline -glob $transfer_from_outside *\;$external_speciation_edge*]
			set transfer_from_outside [lremove $transfer_from_outside $xfer_to_rm]
		} elseif {[lsearch -glob $transfer_from_anoxygeo *\;$edge\)*] > -1} {
			puts "This edge $external_speciation_edge was transferred horizontally from within Anoxy/Geobacillus"

			set events_from_outside [lremove $events_from_outside $external_event]
		}
	}

	puts "Leftover transfer events: [llength $transfer_from_outside]"
	puts "Leftover speciation/duplication events: [llength $events_from_outside]"

	## Any leftover transfer events must be associated with a transfer into a single taxon ##

	foreach external_transfer $transfer_from_outside {
		set final_subevent [get_final_subevent $external_transfer]
		transfer_subevent_parse $final_subevent
		puts "External transfer from \($donor_edge\) to Anoxy/Geobacillus \($receiver_edge\)"
		lappend HGT_table_out "HGT: $donor_edge\t$receiver_edge"
		incr scenario_1_hgt
		set scenario 1
	}


	## As we traverse the branches up, each edge can have a number of events associated with it.
	## Assuming the edges can be followed down to extant species (which is always the case in our bottom up analysis), it is true that any directly ancestral edge must then split into two other edges - UNLESS the edge in question was horizontally transferred from another lineage.
	## Due to the ultrametric nature of the species tree - which is required for the model-based assessment - but incorrect due to equivalent branch lengths assigned to the tree (a quartet based method has no way to assign age-inferred branch lengths). 
	## This arbitrary assignment of branch lengths is interpreted by mowgli as the relative ancestral age, so in reality more broadly represented lineages will typically have higher order branching and thus appear "more ancient" to mowgli. 
	## When inferring direction of HGT, mowgli has to account for this relative age discrepancy - and to make some of these gene trees work, the program will actually predict a more ancestral transfer into a lineage, followed by multiple losses (and perhaps transfers -out- to younger lineages). 
	## This is not an issue, since the closest donor organism will be estimated by parametric methods later; however, an HGT event into an ancestral node followed by losses in all branches apart from the present Anoxy/Geobacillus species should ALSO be considered as a transfer into the Anoxy/Geo taxa ##

	# Now for each speciation that resulted in a topmost anoxy_geo node, we want to check whether the sister lineages have been lost ##
	if {[llength $events_from_outside] != 0} {
		foreach external_event $events_from_outside {
			set final_subevent [get_final_subevent $external_event]
			non_trans_subevent_parse $final_subevent
			## Tested child is refers to the node that results in our anoxy_geo clade ##
			set anoxy_geo_ancestral_branch $edge
			puts "\nThe ancestral edge shared by the anoxy_geo clade is: $anoxy_geo_ancestral_branch"
			set tested_child $child

			## Find all the events where the parent of our speciation is a parent ##
			set children [lsearch -all -inline -glob $events_list *\t\($parent,*]
			foreach child $children {
				if {[regexp "Tran" $child] == 1} {
					set children [lremove $children $child]
				}
			}
			
			## If there is a speciation event, that doesn't correspond to an HGT but also doesn't have a sister group, it must then correspond to the root edge of the tree. In biological terms, this is Mowgli assuming that gene originated within the Anoxy/Geobacillus clade ##

			set children [reduce_non_redundant_children $children]

			if {[llength $children] == 1} {
				puts "Mowgli predicts that this gene originated within the Geobacillus clade"
				lappend HGT_table_out "Is root"
				incr scenario_3_hgt
				set scenario 3
				puts [lindex $children 0]
				continue
			}

			foreach child $children {
				## Isolate the speciation event that is not the one that gives rise to the anoxy_geo clade - i.e. find its sister group ##
				if {[regexp "$parent\,$tested_child" $child] == 0} {
					set other_child $child
					set final_subevent [get_final_subevent $other_child]
					## If the other event is a transfer ... ##
					if {[is_subevent_transfer $final_subevent] == TRUE} {
						transfer_subevent_parse $final_subevent
					## Otherwise set the event type as the fate ##
					} else {
						non_trans_subevent_parse $final_subevent
						set fate $event_type
						set sister_edge $edge
						puts "The fate of the sister group to the ancestral edge $sister_edge == $fate"
					}
				}
			}

			set uplevel 0
			

			while {$fate == "Loss"} {

				incr uplevel
				set tested_child $parent

				## Find events corresponding to the ancestral branch ##
				set up_events [lsearch -all -inline -glob $events_list *,$tested_child\)*]
				

				# puts $up_events
				## If there was only one event, this MUST be the speciation leading down to the known edge ##
				if {[llength $up_events] == 1} {

					set up_event [lindex $up_events 0]

					set final_subevent [get_final_subevent $up_event]
					non_trans_subevent_parse $final_subevent
					set tested_child $child
					set ancestral_edge $edge
					puts "Up-level $uplevel: There was only one ancestral event - corresponding to the edge: $ancestral_edge"

					set children [lsearch -all -inline -glob $events_list *\($parent,*]

					foreach child $children {
						## For all the events that do not include the one we came from ... ##
						if {[regexp "$parent\,$tested_child" $child] == 0} {
							## There can be a number of transfers out before a further speciation, extant, or loss; we need to ignore those ##
							set other_child $child
							set final_subevent [get_final_subevent $other_child]

							if {[is_subevent_transfer $final_subevent] == TRUE} {
								transfer_subevent_parse $final_subevent
								## Differentitate between transfers OUT and transfers IN. A transfer IN can occur if there was a loss in that taxon - i.e. this is a transfer with replacement for that taxon ##
								if {[regexp $parent $donor_edge] == 1} {
									# puts "Transfer from sister group to another branch... ignoring..."
								} else {
									puts "A transfer IN - must be accompanied by a loss in the same branch"
								}
							} else {
								non_trans_subevent_parse $final_subevent
								set fate $event_type
								set sister_edge $edge
							}
						}
					}
					puts "-->The fate of the sister group to the ancestral edge $sister_edge == $fate"

				## If there is more than one event, there could have been a transfer of this edge IN (or a transfer OUT?). If there was a transfer of an edge - and the second event is labelled as a speciation, is the speciation not moot? ##
				} elseif {[llength $up_events] > 1} {
					set event_counter 1
					## Find the new ancestral edge ##
					set ancestral_edge [regexp -inline $pata$tested_child$patb [lindex $up_events 0]]
					puts "Up-level $uplevel: There were [llength $up_events] events - corresponding to the edge: $ancestral_edge"

					## Process the events in turn ##
					foreach event $up_events {
						set final_subevent [get_final_subevent $event]

						if {[is_subevent_transfer $final_subevent] == TRUE} {
							transfer_subevent_parse $final_subevent
							set donor_edge_ends [split $donor_edge \,]
							if {[regexp $tested_child $receiver_edge] == 1 && [lsearch $anoxygeo_present_nodes [lindex $donor_edge_ends 0]] == -1 && [lsearch $anoxygeo_present_nodes [lindex $donor_edge_ends 1]] == -1} {
								puts "\t-->Event $event_counter: This edge \($receiver_edge\) was transferred horizontally from an external source"
								
								lappend HGT_table_out "HGT_multiloss: $donor_edge\t$receiver_edge"

								set fate $trans_type
								break
							}
						} else {
							non_trans_subevent_parse $final_subevent
							set fate $event_type
							set sister_edge $edge
						}
					}
				## If there are no further events, we must have hit the root with all sister branches lost - this is a scenario 4 putative HGT ##	
				} else {
					puts "Mowgli predicts root above AG clade with all other sister branches lost"
					lappend HGT_table_out "Is root with multi-loss"
					set root_multi_loss 1
					break
				}
			}

			if {[regexp "Spec" $fate] == 1 || $fate eq "Dup" || $fate eq "Extant"} {
				puts "The sister branch - $sister_edge - further speciated, so the Anoxy/Geobacillus clade subtended by $anoxy_geo_ancestral_branch was NOT transferred"
			} elseif {[regexp "Tran" $fate] == 1} {
				puts "The ancestral edge $ancestral_edge transferred horizontally, and all other clades descendants apart from the Anoxy/Geobacillus clade (subtended by $anoxy_geo_ancestral_branch) were lost. This is a putative transfer."
				incr scenario_2_hgt
				set scenario 2
			} elseif {$root_multi_loss == 1} {
				incr scenario_4_hgt
				set scenario 4
			} else {
				puts "Hmmmm"
				exit 2
			}

		}
	}

	if {[llength $HGT_table_out] == 1} {
		lappend HGT_table_out "No transfers into Anoxy/Geobacillus"
		lappend global_hgt_table "$out_dir\t0\t$scenario"
		incr vertical
	} else {
		lappend global_hgt_table "$out_dir\t1\t$scenario"
		incr horizontal
	}
	# Write out the file containing all the donor and receiver edges for the transfers INTO geobacillus from an external source ##
	set out [open Transfers_in.tsv w]
	puts $out [join $HGT_table_out \n]
	close $out

	puts stdout [join $HGT_table_out \n]

	incr all_tested
}

set global_hgt_table [lsort -dictionary $global_hgt_table]

if {$ignore_inside_only_switch == 0} {
	set out_table_filename $direct/Results_Full/T$transfer_cost\_D2_L1_results.txt
	set out_stats_filename $direct/Results_Full/T$transfer_cost\_D2_L1_stats.txt
} else {
	set out_table_filename $direct/Results_NO_IG/T$transfer_cost\_D2_L1_NO_IG_results.txt
	set out_stats_filename $direct/Results_NO_IG/T$transfer_cost\_D2_L1_NO_IG_stats.txt
}

set global_hgt_table [join $global_hgt_table \n]
set out [open $out_table_filename w]
puts $out $global_hgt_table
close $out

set stats [open $out_stats_filename w]
puts $stats "Total number of gene families tested: $all_tested\nAG vertical: $vertical\nAG horizontal: $horizontal\n\nScenario 1: $scenario_1_hgt\nScenario 2: $scenario_2_hgt\nScenario 3: $scenario_3_hgt\nScenario 4: $scenario_4_hgt"
close $stats