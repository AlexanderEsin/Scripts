#!/usr/bin/tclsh

### This script splits the accumulated blast hit files (in out_all) into
### a number of folders spcecified by n. The blast hits are paired 
### for downstream RBBH analysis.

## Proc to find the opposite pair for a blast out file
proc RearrangeBlastName {blast_file_name} {
	set split_name	[split [string range $blast_file_name 0 end-4] \&]
	set reassemble	"[lindex $split_name 1]&[lindex $split_name 0]"
	set with_suffix	"$reassemble\.tsv"
	return $with_suffix
}

set group_name	[lindex $argv 0]

## Set the directories and make the output_split master dir
set direct		/scratch/ade110/Geo_again/Consensus_groups/$group_name/BLAST
set blast_out	$direct/OUT_all_split
set out_split	$direct/OUT_all_rearrange
set temp_hold	$direct/temp
file mkdir $out_split $temp_hold

## Copy all the blast outputs into a global temp folder
set input_out_dirs	[glob -type d $blast_out/*]
foreach input_dir $input_out_dirs {
	set blast_out_files	[glob $blast_out/$input_dir/*tsv]
	foreach blast_file $blast_out_files {
		file copy -force $blast_file $temp_hold
	}
}

## Check how many files are now in the temp dir
set temp_files	[glob $temp_hold/*tsv]
puts stdout		"We moved [llength $temp_files] BLAST output files into a temporary holding director $temp_hold"


## // Now split these files into subdirectories according to name // ##

## Make the correct number of folders for the split
# n = number of split folders wanted
set n 64
set i 1
while {$i <= $n} {
	file mkdir $out_split/$i
	incr i
}

## Get a list of all the blast output files
cd $temp_hold
set blast_out_files [glob *tsv]
set num_out_files	[llength $blast_out_files]

## Set up the variables to keep track of the splitting
set folder_counter	0
set pair_counter	0
set files_remain	$num_out_files

## While files remain, sort each file and its RBBH pair into
## split folders
while {$files_remain > 0} {
	
	# If we've put a pair of blast results into each folder
	# Reset and start filling at the start
	if {$folder_counter == $n} {
		set folder_counter 0
	}

	# Get the right directory
	set split_dir	[expr $folder_counter + 1]
	# Pick a file and get its RBBH opposite
	set file_A		[lindex $blast_out_files $pair_counter]
	set file_B		[RearrangeBlastName $file_A]

	# If we can't find file_A, then it could have already
	# been copied as file_B
	if {[file exists $file_A] != 1} {
		incr pair_counter
		continue
	}
	# If we can't find file_B this must be an error
	if {[file exists $file_B] != 1} {
		error "Can't find $file_B"
	}

	# Move the files
	file rename $file_A	$out_split/$split_dir
	file rename $file_B	$out_split/$split_dir
	
	# Incr the pair and folder counter
	incr pair_counter
	incr folder_counter
	# Counter to stop at the right time
	set files_remain	[expr $files_remain - 2]
}

cd $direct
file delete -force $temp_hold