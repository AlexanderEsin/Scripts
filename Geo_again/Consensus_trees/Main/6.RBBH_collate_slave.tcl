#!/usr/bin/tclsh

source /home/ade110/Scripts/General_utils.tcl

set group_name	[lindex $argv 0]
set eval_cutoff	[lindex $argv 1]
set eval_num	[string range $eval_cutoff 2 end]

set direct			/scratch/ade110/Geo_again/Consensus_groups/$group_name/RBBH
set RBBH_dir		$direct/RBBH_$eval_num

## The inparalog directory should exist from Find_true_paralogs.tcl
set para_RBBH_dir	$direct/InParalogs_RBBH

## Define ortholog and master RBBH dirs: then make them
set orth_RBBH_dir	$direct/Orthologs_RBBH
set master_RBBH_dir	$direct/Master_RBBH
file mkdir $orth_RBBH_dir $master_RBBH_dir


## If the ortholog (because we append) file already exists: delete
if {[file exists $orth_RBBH_dir/Orthologs_RBBH_$eval_num\.txt] == 1} {
	file delete	$orth_RBBH_dir/Orthologs_RBBH_$eval_num\.txt
}

## Get a list of all the RBBH files
cd $RBBH_dir
set RBBH_file_list [glob *.txt]
set num_RBBH_files	[llength $RBBH_file_list]

## Prepare a log file to track progress
set logchan			[open $master_RBBH_dir/Master_progress_$eval_num\.log a]

## Add the data body of each RBBH file to an Ortholog file
set counter 1
foreach RBBH_file $RBBH_file_list {

	set RBBH_data	[string trim [openfile $RBBH_file]]
	set no_header	[string trim [string range $RBBH_data [string first \n $RBBH_data] end]]

	set out			[open $orth_RBBH_dir/Orthologs_RBBH_$eval_num\.txt a]
	puts $out		$no_header
	close $out

	puts $logchan	"Added $counter / $num_RBBH_files RBBH files to Orthologs_RBBH"
	incr counter
}


puts $logchan	"Adding Ortholog and InParalog RBBHs into Master ..."
## Add the ortholog and inParalog RBBHs into master
if {[file exists $para_RBBH_dir/InParalogs_RBBH_$eval_num\.txt] == 1} {

	# Set our out file
	set out		[open $master_RBBH_dir/Master_RBBH_$eval_num\.txt w]
	# Copy in orthologs
	set in		[open $orth_RBBH_dir/Orthologs_RBBH_$eval_num\.txt]
	fcopy		$in $out
	close $in
	# Copy in paralogs
	set in		[open $para_RBBH_dir/InParalogs_RBBH_$eval_num\.txt]
	fcopy		$in $out
	close $in
	# Close our out file
	close $out

} else {
	error "Cannot find correct InParalog RBBH file ..."
}

puts $logchan	"Master RBBH file created"
close $logchan