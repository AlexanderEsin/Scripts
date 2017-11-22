#!/usr/bin/tclsh

## Note the shebang for AX4 above
## In this script we run FastTree on each of the selected
## alignments in parallel. The FastTree/Split_MSA folder was
## tar-balled and rsynced to the ax4 Scratch into the
## $SCRATCH/FastTree directory. There it was untarred

source		/home/ade110/Scripts/General_utils.tcl

set direct			/scratch/ade110/FastTree/Split_MSA
set all_tree_dir	$direct/../All_trees
set split			[lindex $argv 0]
# set split	1

set align_folders	[glob -type d $direct/$split/*]

foreach align_folder $align_folders {
	set group_number	[file tail $align_folder]

	set input_align		$align_folder/$group_number\.phy
	set output_tree		$align_folder/FT_$group_number\_tree.txt

	puts stdout			"Building FT tree for alignment: $group_number in $align_folder"

	## If the tree file doesn't already exist, build tree
	## with FastTree
	if {[file exists $output_tree] != 1} {
		set FT_log			[open $align_folder/FT_log.txt a]
		catch {exec FastTree -out $output_tree $input_align >&@$FT_log}
		close $FT_log
	}

	## We want to relabel the tips of the tree
	## (Output from fasttree is in k[key] format)
	set key_file		$align_folder/KEY_$group_number\.txt
	set key_data		[split [string trim [openfile $key_file]] \n]

	set tree_data		[string trim [openfile $output_tree]]
	set tree_relab		$tree_data

	set label_pat		{k+?[0-9]+?:+?}
	set tree_leaves		[regexp -all -inline $label_pat $tree_data]

	## For each tree tip, rename it according to old nomenclature
	## for compliance downstream.
	foreach tip $tree_leaves {
		# Trim away lagging colon
		set tip_trim	[string range $tip 0 end-1]
		set key_entry	[split [lsearch -inline $key_data "$tip_trim\t*"] \t]

		set binomial	[lindex $key_entry 1]
		set class		[lindex $key_entry 2]
		set id			[lindex $key_entry 3]

		set new_tip		"$binomial \{$class\} \{$id\}:"
		set tree_relab	[string map [list $tip $new_tip] $tree_relab]	
	}

	set out		[open $align_folder/FT_$group_number\_relab.txt w]
	puts $out	$tree_relab
	close $out

	## Finally copy each relabelled tree into a joint directory
	## for easy extraction.
	file copy -force $align_folder/FT_$group_number\_relab.txt $all_tree_dir/FT_$group_number\_relab.txt
}

puts stdout				"\n\nBuilding FT trees for Split_MSA folder: $split complete"