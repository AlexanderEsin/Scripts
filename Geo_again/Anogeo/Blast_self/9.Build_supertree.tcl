#!/usr/local/bin/tclsh

## Here we build a RAxML consensus tree using the
## supermatrix alignment. We allow RAxML to determine
## the best protein model automatically.
## Using: RAxML v8.2.10
source ~/Documents/Scripts/General_utils.tcl

set threads			20
# set evalue			1e-150
set evalue			[lindex $argv 0]
set trunc_eval		[string range $evalue 2 end]

## Define the directories and create new output folder
set direct			/users/aesin/Desktop/Geo_again/Anogeo_analysis

set groups_dir		$direct/Family_groups/$trunc_eval
set super_phy_file	$groups_dir/Supermatrix_phylip.phy

set output_dir		$groups_dir/RAxML_tree
file mkdir			$output_dir

## Check the output folder is empty, as RAxML throws errors
if {[llength [glob $output_dir/*]] > 0} {
	set tree_files		[glob -type f $output_dir/*]
	foreach file $tree_files {
		file delete $file
	}
	puts "\nOutput folder cleaned of [llength $tree_files] files!\n"
}

## // Build the RAxML tree // ##
## The "-" of the evalue is not allowed in the run ID
## so trimming to just the exponential
set output_name 	"super_tree$trunc_eval.txt"
set random_pseed	[expr int(100000 * rand())]
set random_xseed	[expr int(100000 * rand())]
set bootstraps		100

## Record the options in a log file - including seeds
set out		[open $output_dir/parameters.txt w]
puts $out	[join [concat "evalue=$trunc_eval" "input=$super_phy_file" "out_name=$output_name" "p=$random_pseed" "x=$random_xseed" "BS=$bootstraps" "threads=$threads"] \n]
close $out

## Build tree
puts "Building RAxML super tree at e-value: $trunc_eval ..."
catch {exec raxml -f a -s $super_phy_file -n $output_name -w $output_dir -m PROTCATAUTO -p $random_pseed -x $random_xseed -N $bootstraps -T $threads}

