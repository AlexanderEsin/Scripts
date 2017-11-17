#!/usr/local/bin/tclsh

source ~/Documents/Scripts/General_utils.tcl

set eval_cutoff	[lindex $argv 0]
set eval_num	[string range $eval_cutoff 2 end]

set direct			/users/aesin/Desktop/Geo_again/Anogeo_analysis/RBBH
set RBBH_dir		$direct/RBBH_$eval_num
## The inparalog directory should exist from 2.Find_true_paralogs.tcl
set para_RBBH_dir	$direct/InParalogs_RBBH

## Define ortholog and master RBBH dirs: then make them
set orth_RBBH_dir	$direct/Orthologs_RBBH
set master_RBBH_dir	$direct/Master_RBBH
file mkdir $orth_RBBH_dir $master_RBBH_dir


## If the master file already exists: delete
if {[file exists $orth_RBBH_dir/Orthologs_RBBH_$eval_num\.txt] == 1} {
	file delete	$orth_RBBH_dir/Orthologs_RBBH_$eval_num\.txt
}

## Get a list of all the RBBH files
cd $RBBH_dir
set RBBH_file_list [glob *.txt]

## Add the data body of each RBBH file to an Ortholog file
set counter 1
foreach RBBH_file $RBBH_file_list {

	set RBBH_data	[string trim [openfile $RBBH_file]]
	set no_header	[string trim [string range $RBBH_data [string first \n $RBBH_data] end]]

	set out			[open $orth_RBBH_dir/Orthologs_RBBH_$eval_num\.txt a]
	puts $out		$no_header
	close $out

	puts "Added $counter/[llength $RBBH_file_list] of RBBH files to Orthologs_RBBH"
	incr counter
}


puts "Adding Ortholog and InParalog RBBHs into Master ..."
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

puts "Master RBBH file created"