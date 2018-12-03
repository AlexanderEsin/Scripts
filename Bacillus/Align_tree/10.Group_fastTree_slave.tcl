source /home/ade110/Scripts/General_utils.tcl

set job_number		[lindex $argv 0]
set index_number	[expr $job_number - 1]

puts $job_number

## Directories
set direct 			/home/ade110/Work/Bacillus
set align_fin_dir	$direct/Group_alignment/Final_aligned
set tree_bin_dir	$direct/Group_fastTree/Bin
set tree_fin_dir	$direct/Group_fastTree/Final_trees
file mkdir			$tree_bin_dir $tree_fin_dir


## A list of the alignments, select the one for
## this subjob
set final_aligns	[glob $align_fin_dir/*fas]
set group_align		[lindex $final_aligns $index_number]

## Create a seperate directory for this alignment
## Reconstruct tree within that directory (bin) so we can output log
## Final trees are then put into the final tree directory
set file_name	[file tail $group_align]
set group_num	[string range $file_name 0 [string first \_ $file_name]-1]

set group_out_dir	$tree_bin_dir/$group_num
file mkdir			$group_out_dir

## Set up + open log channel
set log_file		"$group_num\_FT_log.txt"
if {[file exists $group_out_dir/$log_file] == 1} {
	file delete $group_out_dir/$log_file
}
set FT_log			[open $group_out_dir/$log_file a]

## Run FastTree - gamma correction at end and use standard bootstraps
## (instead of default SH test). Using standard JTT model for AA
catch {exec FastTree -gamma -boot 100 -out $group_out_dir/$group_num\_FT_tree.txt $group_align >&@$FT_log}

## Close log channel
close $FT_log

## If tree is written out, copy to final directory
if {[file exists $group_out_dir/$group_num\_FT_tree.txt] == 1} {
	file copy -force $group_out_dir/$group_num\_FT_tree.txt $tree_fin_dir/$group_num\_FT_tree.txt
}