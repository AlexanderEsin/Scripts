#!/usr/local/bin/tclsh

## In this script we split the alignment+key containing
## directories into sub-sections to be assembled in
## parallel on the AX4.

set direct		/Users/aesin/Desktop/FastTree

## Get the input directories containing keys + alignments
set input_dir	$direct/Input_groups
set group_dirs	[glob -type d $input_dir/*]
set num_dirs	[llength $group_dirs]

## Make the output directories
set out_split	$direct/Split_MSA
file mkdir		$out_split

## Make the correct number of folders for the split
# n = number of split folders wanted
set n 50
set i 1
while {$i <= $n} {
	file mkdir $out_split/$i
	incr i
}


## Set up the variables to keep track of the splitting
set align_counter	0
set folder_counter	0
set dirs_remain		$num_dirs

## While files remain, sort each file and its RBBH pair into
## split folders
while {$dirs_remain > 0} {
	
	# If we've put a directory into each folder,
	# reset and start filling at the start
	if {$folder_counter == $n} {
		set folder_counter 0
	}

	# Get the right directory
	set split_dir	[expr $folder_counter + 1]

	# Pick the next directory
	set group_dir	[lindex $group_dirs $align_counter]

	# Move the files
	file copy -force $group_dir	$out_split/$split_dir
	
	# Incr the pair and folder counter
	incr align_counter
	incr folder_counter
	# Counter to stop at the right time
	set dirs_remain	[expr $dirs_remain - 1]
}

