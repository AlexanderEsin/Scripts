#!/usr/local/bin/tclsh
## Run mowgli on the HPC ##

## For desktop ##
source 		~/Documents/Scripts/General_utils.tcl

set direct		/users/aesin/Desktop/FastTree
set input_dir	$direct/All_trees_relab
set out_all 	$direct/Mowgli_input


set penalties [list 3 4 5 6]


foreach penalty $penalties {

	set out_penalty $out_all/Input_$penalty
	file mkdir $out_penalty

	## Get list of all tree files in the input directory
	set trees [lsort -dictionary [glob $input_dir/*txt]]

	set folder_counter 1
	## Foreach tree, make a new directory to house the tree file
	foreach tree $trees {
		file mkdir	$out_penalty/$folder_counter
		file copy -force $tree $out_penalty/$folder_counter

		incr folder_counter
	}

}



