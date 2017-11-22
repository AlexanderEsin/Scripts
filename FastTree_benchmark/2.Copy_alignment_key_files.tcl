#!/usr/local/bin/tclsh

## In this script we transfer the alignments and key files
## for the selected groups into a the testing directory. We will
## split these files into folders for FastTree reconstruction
## on the AX4

source ~/Documents/Scripts/General_utils.tcl

set direct			/Users/aesin/Desktop
set align_files		$direct/Geo_analysis/Geo_v_all/2.0/Tree_build_final

set output_dir		$direct/FastTree/Input_groups
set group_file_list	$output_dir/Group_list.txt


## Get a list of groups for which we need alignments + keys
set groups_needed	[split [string trim [openfile $group_file_list]] \n]


## For each group copy the alignment and key into the FastTree testing directory
foreach group $groups_needed {
	# Make output directory
	set group_output	$output_dir/$group
	file mkdir			$group_output

	# Define alignment and key files
	set alignment_file	"$align_files/$group/$group\.fas"
	set key_file		"$align_files/$group/KEY_$group\.txt"

	# Copy them
	file copy -force $alignment_file $key_file $group_output

	## Remake the phylip format. Some of the phylip
	## alignments appear to be corrupted (e.g. indentation in 1432)
	## Seqret reads them fine but FastTree throws an error

	catch {exec seqret $group_output/$group\.fas $group_output/$group\.phy -osformat phylip 2> /dev/null}
}

puts "Copying alignment and key files ... DONE"