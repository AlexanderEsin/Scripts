#!/usr/local/bin/tclsh
## Rename tree labels from the key ##

source "/Users/aesin/Documents/Scripts/General_utils.tcl"

## Find and rein in tree and key files ##
set tree_file "/Users/aesin/Desktop/Geo_analysis/New_geobacillus_consensus/Tree/RAxML_bipartitions.concat_geo_tree.txt"
set key_file "/Users/aesin/Desktop/Geo_analysis/New_geobacillus_consensus/Master_alignment_key.txt"

set tree_data	[string trim [openfile $tree_file]]
set key_data	[split [string trim [openfile $key_file]] \n]

## For each entry in the key, find and rename the corresponding tip label in the tree ##

foreach key_entry $key_data {
	set key_entry_l		[split $key_entry \t]
	set key_value		[lindex $key_entry_l 0]
	set key_binomial	[lindex $key_entry_l 1]

	regsub -all "$key_value\:" $tree_data "$key_binomial\:" tree_data
}


set out [open "/Users/aesin/Desktop/Geo_analysis/New_geobacillus_consensus/Tree/Final_consensus_out.txt" w]
puts $out [string trim $tree_data]
close $out
