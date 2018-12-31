# source /home/ade110/Scripts/General_utils.tcl
source /mnt/storage/home/aesin/Scripts/General_utils.tcl

set job_number		[lindex $argv 0]
set index_number	[expr $job_number - 1]

## Directories
# set direct 			/home/ade110/Work/Staph
set direct 			/mnt/storage/rawdata/Warnecke_Rawdata/ADE/Work/Staph
set red_grp_dir		$direct/Family_mapping/Reduced_groups_fasta
set align_bin_dir	$direct/Group_alignment/Bin
set align_fin_dir	$direct/Group_alignment/Final_aligned
file mkdir			$align_bin_dir $align_fin_dir

## A list of all the reduced fasta files
set final_group_fastas	[glob $red_grp_dir/*fasta]
set group_fasta_file	[lindex $final_group_fastas $index_number]

## For the group fasta file, create a seperate directory
## Align within that directory (bin) so we can output log
## Final alignments are then put into the final align dir
set file_name	[file tail $group_fasta_file]
set group_num	[string range $file_name 0 end-6]

set group_out_dir	$align_bin_dir/$group_num
file mkdir			$group_out_dir

## If the output alignment already exists in the final directory, exit
if {[file exists $align_fin_dir/$group_num\_aligned.fas] == 1} {
	exit
}

## Otherwise run clustalo alignment and write the final alignment to
## the final directory
catch {exec clustalo --auto -i $group_fasta_file -o $group_out_dir/$group_num\_aligned.fas --outfmt=fasta -v --force -l $group_out_dir/$group_num\_MSA_log.txt --threads=1}
if {[file exists $group_out_dir/$group_num\_aligned.fas] == 1} {
	file copy -force $group_out_dir/$group_num\_aligned.fas $align_fin_dir/$group_num\_aligned.fas
}
