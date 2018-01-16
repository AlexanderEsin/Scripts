#!/usr/local/bin/tclsh

## Here we build a FastTree consensus tree using the
## supermatrix alignment. We make FastTree use the LG
## protein model automatically (RAxML estimated)
## Using: FastTree 2.1.10
source ~/Documents/Scripts/General_utils.tcl

set threads			8
set group_name		[lindex $argv 0]
set evalue			[lindex $argv 1]
set trunc_eval		[string range $evalue 2 end]

## Define the directories and create new output folder
set direct			/users/aesin/Desktop/Geo_again
set consens_dir		$direct/Consensus_groups/$group_name

set groups_dir		$consens_dir/Family_groups/$trunc_eval
set super_phy_file	$groups_dir/Supermatrix_phylip.phy

set output_dir		$groups_dir/FT_tree
file mkdir			$output_dir

## Check the output folder is empty; if not clear
if {[llength [glob -nocomplain $output_dir/*]] > 0} {
	set tree_files		[glob -type f $output_dir/*]
	foreach file $tree_files {
		file delete $file
	}
	puts "\nOutput folder cleaned of [llength $tree_files] files!\n"
}

## // Build the FastTree tree // ##
## The "-" of the evalue is not allowed in the run ID
## so trimming to just the exponential
set output_name 	"FT_super_tree$trunc_eval.txt"

## Record the options in a log file - including seeds
set FT_log		[open $output_dir/FT_log.txt w]

## Build tree
puts "Building FastTree super tree at e-value: $trunc_eval ..."
catch {exec FastTree -lg -out $output_dir/$output_name $super_phy_file >&@$FT_log}
close $FT_log