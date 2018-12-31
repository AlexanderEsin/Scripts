#!/usr/local/bin/tclsh
source ~/Documents/Scripts/General_utils.tcl
package require sqlite3

## Paths and output folders
set master_dir		/Users/aesin/Desktop/Bacillus

# Find and open the DB
set all_db_file		$master_dir/ForBac_prot_db
sqlite3 allProt_db	$all_db_file

# Directory to find tipID <-> protID keys
set tipKey_dir		$master_dir/Mowgli/GeneTree_input

# Directories to read HGT and vertical events
set mowgli_dir		$master_dir/Mowgli/Mowgli_output
set refined_dir		$mowgli_dir/Refined_events/Events
set vertical_dir	$mowgli_dir/Constant_events

# Define and make output directory
set cleaned_dir		$mowgli_dir/Cleaned_events
set hgtClean_dir	$cleaned_dir/HGT_events
set verClean_dir	$cleaned_dir/Ver_events
file mkdir			$hgtClean_dir $verClean_dir

## Penalty and HGT-type lists
# set penalty_list	[list 3 4 5 6]
set penalty_list	[list 3 4]
set hgtType_list	[list "lHGT" "sHGT" "Ver"]

## Open log file to record cleaning results
set cleanLog_chan	[open $cleaned_dir/HGT_cleaningLog.log w]


## Core taxids
set core_taxid_file		$master_dir/Core_genomes/Genome_lists/coreToKeep.tsv
set core_taxid_tbl		[lrange [split [string trim [openfile $core_taxid_file]] \n] 1 end]
set core_taxids		{}
foreach row $core_taxid_tbl {
	set taxid	[lindex [split $row \t] 2]
	lappend core_taxids $taxid
}


## Iterate over both types of HGT and all penalties
foreach hgtType $hgtType_list {

	foreach penalty $penalty_list {

		if {$hgtType ne "Ver"} {
			set input_file	$refined_dir/T$penalty\_full_$hgtType\_events.tsv
			set refined_data	[split [string trim [openfile $input_file]] \n]

			set refinedHeader	[lindex $refined_data 0]
			set refined_data	[lrange $refined_data 1 end]
		} else {
			set input_file		[glob $vertical_dir/$hgtType\_const_t$penalty*]
			set refined_data	[split [string trim [openfile $input_file]] \n]
		}

		# Counters
		set process_counter		1
		set mixedEvent_counter	0
		set cleanEvent_counter	0
		set plasEvent_counter	0

		# Output lists
		set mixedEvent_list		{}
		set cleanEvent_list		{}
		set plasEvent_list		{}

		# Process each refined event
		foreach event $refined_data {
			# Keep track of process
			puts -nonewline	"\rProcessing refined event $process_counter // [llength $refined_data] ..."
			flush stdout

			if {$hgtType ne "Ver"} {
				set columns		[split $event \t]
				set groupNum	[lindex $columns 1]
				set tipID_list	[split [lindex $columns end] " "]
				
				# Open the key to translate tipID to protID
				set key_file	$tipKey_dir/$groupNum/$groupNum\_KEY_tips.txt
				set key_data	[lrange [split [string trim [openfile $key_file]] \n] 1 end]; # Remove header as well

				# For each tip in the event, get the protID and check whether it's plasmid-borne
				set isPlasmid_list	{}
				foreach tipID $tipID_list {
					set key_entry	[lsearch -all -glob -inline $key_data "*\t$tipID"]
					if {[llength $key_entry] == 1} {set key_entry [lindex $key_entry 0]} else {error "Finding too many key entries for $tipID in $groupNum at $penalty at $hgtType"}
					set protID		[lindex [split $key_entry \t] 0]
					
					set isPlasmid	[allProt_db eval {SELECT plasmid FROM t1 WHERE protID = $protID}]
					lappend isPlasmid_list $isPlasmid
				}

			} else {
				set isPlasmid_list {}
				foreach core_taxid $core_taxids {
					set isPlasmid	[allProt_db eval {SELECT plasmid FROM t1 WHERE OrthGroup = $event AND taxid = $core_taxid}]
					if {[string length $isPlasmid] != 0} {
						set unique [lsort -unique $isPlasmid]
						lappend isPlasmid_list $unique
					}
				}
			}
			

			# Check we have no NAs (all AG protIDs should have a T or F for plasmid)
			if {[llength [lsearch -all $isPlasmid_list "NA"]] > 0} {error "Finding NAs for isPlasmid for $event at $penalty at $hgtType"}

			# Get just unique plasmid states (all T || all F || T + F)
			set isPlasmid_uniq	[lsort -unique $isPlasmid_list]

			if {[llength $isPlasmid_uniq] == 2} {
				incr mixedEvent_counter
				lappend mixedEvent_list	$event
			} elseif {[llength $isPlasmid_uniq] > 2 || [llength $isPlasmid_uniq] == 0} {
				error "Wrong number of isPlasmid fates \($isPlasmid_uniq\) for $event at $penalty at $hgtType"
			} elseif {$isPlasmid_uniq eq "F"} {
				incr cleanEvent_counter
				lappend cleanEvent_list	$event
			} elseif {$isPlasmid_uniq eq "T"} {
				incr plasEvent_counter
				lappend plasEvent_list	$event
			} else {
				error "Something wrong with \($isPlasmid_uniq\) for $event at $penalty at $hgtType"
			}

			incr process_counter
		}

		# Prepare the clean events to be written out
		if {$hgtType ne "Ver"} {
			set cleanEvent_data	[concat [list [join $refinedHeader \t]] $cleanEvent_list]
			set cleanOut_dir	$hgtClean_dir
		} else {
			set cleanEvent_data $cleanEvent_list
			set cleanOut_dir	$verClean_dir
		}

		set cleaned_file	$cleanOut_dir/[file tail $input_file]
		
		# Write out clean events
		set cleaned_out		[open $cleaned_file w]
		puts $cleaned_out	[join $cleanEvent_data \n]
		close $cleaned_out

		# Write stats to log file
		multiputs $cleanLog_chan stdout "\nResults for $hgtType at penalty = $penalty:"
		multiputs $cleanLog_chan stdout "Total events processed: [llength $refined_data] of which ..."
		multiputs $cleanLog_chan stdout "\t$cleanEvent_counter events are plasmid-free"
		multiputs $cleanLog_chan stdout "\t$plasEvent_counter events are all plasmid-borne"
		multiputs $cleanLog_chan stdout "\t$mixedEvent_counter events are mixed between chromosome and plasmid\n\n"

	}
}

## Close log file
close $cleanLog_chan

## Close database
allProt_db close
