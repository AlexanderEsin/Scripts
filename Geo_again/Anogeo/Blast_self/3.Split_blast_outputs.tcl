#!/usr/local/bin/tclsh

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

## Set the directories and make the output_split master dir
set direct		/users/aesin/Desktop/Geo_again/Anogeo_analysis
set blast_out	$direct/Blast_anogeo/Blast_out
set out_split	$direct/Blast_anogeo/Blast_out_split
file mkdir $out_split

## Make the correct number of folders for the split
# n = number of split folders wanted
set n 20
set i 1
while {$i <= $n} {
	file mkdir $out_split/$i
	incr i
}

## Get a list of all the blast output files
cd $blast_out
set blast_out_files [glob *.tsv]
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



########## Refill incase of mistake ##########

# cd /scratch/ade110/Consensus_trees/Clostridia/Out_split
# set i 1
# set dirs [glob -type d *]
# foreach dir $dirs {
# 	cd /scratch/ade110/Consensus_trees/Clostridia/Out_split/$dir
# 	set files [glob -nocomplain *.tsv]
# 	foreach file $files {
# 		file rename $file /scratch/ade110/Consensus_trees/Clostridia/Out_all
# 		puts "$i"
# 		incr i
# 	}
# }

# set total_length 0
# set dirs [glob -type d *]
# foreach dir $dirs {
# 	cd /scratch/ade110/Consensus_trees/Bacillaceae/Out_split/$dir
# 	set files [glob -nocomplain *.tsv]
# 	set total_length [expr $total_length + [llength $files]]
# }
# puts "NO: $total_length"

# cd /csc/results/Warnecke/ADE/Backup/Consensus_trees/Bacillales/Out_all
# set files_1 [glob *tsv]
# set out_all_length [llength $files_1]

# cd /scratch/ade110/Consensus_trees/Bacillales/Out_split
# set i 1
# set dirs [glob -type d *]
# set total_length 0
# set missing {}
# foreach dir $dirs {
# 	cd /scratch/ade110/Consensus_trees/Bacillales/Out_split/$dir
# 	set files [glob -nocomplain *.tsv]
# 	set total_length [expr $total_length + [llength $files]]
# 	foreach file $files {
# 		if {[file exists /csc/results/Warnecke/ADE/Backup/Consensus_trees/Bacillales/Out_all/$file] != 1} {
# 			lappend missing "$dir === $file"
# 			break
# 		}
# 		puts "$i"
# 		incr i
# 	}
# }

# puts "\n\n$missing"
# puts $total_length
# puts $out_all_length

