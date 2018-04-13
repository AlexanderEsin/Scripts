#!/usr/local/bin/tclsh

source ~/Documents/Scripts/General_utils.tcl
package require sqlite3


# Paths and output folders
set master_dir		/Users/aesin/Desktop/Geo_again
set baseGroup_dir	$master_dir/Family_groups/-10
set hgtPosition_dir	$master_dir/HGT_position

# Open the basal final groups file
set finalGroup_file	$baseGroup_dir/Final_groups.tsv
set finalGroup_list	[split [string trim [openfile $finalGroup_file]] \n]

# The dnaA protein group is "140"
set dnaAGroup_entry	[lsearch -glob -inline $finalGroup_list "140\t*"]
set dnaAgroup_prots	[lrange [split $dnaAGroup_entry \t] 1 end]

set out		[open $hgtPosition_dir/dnaA_protIDs.txt w]
puts $out	[join $dnaAgroup_prots \n]
close $out
