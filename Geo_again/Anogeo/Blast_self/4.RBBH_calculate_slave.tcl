#!/usr/local/bin/tclsh

source ~/Documents/Scripts/General_utils.tcl

## Proc to find the opposite pair for a blast out file
proc RearrangeBlastName {blast_file_name} {
	set split_name	[split [string range $blast_file_name 0 end-4] \&]
	set reassemble	"[lindex $split_name 1]&[lindex $split_name 0]"
	set with_suffix	"$reassemble\.tsv"
	return $with_suffix
}

## The script argument - a number - defines the input blast folder
set eval_cutoff		[lindex $argv 0]
set folder_number	[lindex $argv 1]

## Set the input and output directories
set direct		/users/aesin/Desktop/Geo_again/Anogeo_analysis
set blast_split	$direct/Blast_anogeo/Blast_out_split

set RBBH_dir	$direct/RBBH

## A cut off value for which results to pick from the blast output
## the default blast has a cut off of 1e-10
set eval_num	[string range $eval_cutoff 2 end]

set RBBH_out	$RBBH_dir/RBBH_$eval_num
set para_out	$RBBH_dir/Paralog_comparison
file mkdir		$RBBH_out; file mkdir $para_out


## Get a list of the input files
cd $blast_split/$folder_number
set blast_file_list	[glob *.tsv]
set num_blast_files	[llength $blast_file_list]


set files_remaining $num_blast_files
while {$files_remaining > 0} {

	set file_A	[lindex $blast_file_list 0]
	set file_B	[RearrangeBlastName $file_A]

	# Check the reciprocal blast output is here
	if {[file exists $file_B] != 1} {
		error "Cannot find $file_B"
	}

	set rbbh_output			{"GeneA\tGeneB\tEval"}
	set paralog_comparison	{"GeneA\tGeneB\tBits"}
	set rbbh_counter 0
	set below_cutoff 0
	set miss_counter 0
	set hits_processed 0

	set hits_A	[split [string trim [openfile $file_A]] \n]
	set hits_B	[split [string trim [openfile $file_B]] \n]

	foreach blast_hit $hits_A {
		set cols_A	[split $blast_hit \t]
		set gene_1	[lindex $cols_A 0]
		set gene_2	[lindex $cols_A 1]

		set RBBH_hit	[lsearch -all -inline $hits_B "$gene_2\t$gene_1*"]

		if {[llength $RBBH_hit] == 1} {
			set hit		[lindex $RBBH_hit 0]

			set eval_A	[lindex $cols_A 10]
			set eval_B	[lindex [split $hit \t] 10]

			set bits_A	[lindex $cols_A 11]
			set bits_B	[lindex [split $hit \t] 11]
			set av_bits	[expr ($bits_A + $bits_B) / 2]

			## A condensed format for in-paralog comparison
			## using bitscore
			lappend paralog_comparison	"$gene_1\t$gene_2\t$av_bits"

			## If the evalue cutoff is not at default minimum for blast (1e-10)
			## check that it's below the treshold both ways
			if {$eval_cutoff == 1e-10} {
				lappend rbbh_output "$gene_1\t$gene_2\t$eval_A"
				lappend rbbh_output "$gene_2\t$gene_1\t$eval_B"
				incr rbbh_counter
			} else {
				## If it's not below threshold, count it separately
				if {$eval_A < $eval_cutoff && $eval_B < $eval_cutoff} {
					lappend rbbh_output "$gene_1\t$gene_2\t$eval_A"
					lappend rbbh_output "$gene_2\t$gene_1\t$eval_B"
					incr rbbh_counter
				} else {
					incr below_cutoff
				}
			}
		## If the Blast hit has no reciprocal...
		} elseif {[llength $RBBH_hit] == 0} {
			incr miss_counter
		## If the blast hit has multiple reciprocals...
		} else {
			error "Don't expect to have more than one reciprocal hit with current settings:\n$blast_hit\n$RBBH_hit"
		}

		incr hits_processed
	}

	## Preapre the output file name for the RBBH hits
	set RBBH_out_name	[string range $file_A 0 end-4]

	## Write out the RBBH hits
	set out		[open $RBBH_out/$RBBH_out_name\.txt w]
	puts $out	[join $rbbh_output \n]
	close $out

	## Write out the condensed RBBH with bitscores for comparison
	## to in-paralogs
	set out		[open $para_out/$RBBH_out_name\.txt w]
	puts $out	[join $paralog_comparison \n]
	close $out

	## Counter for loop
	set files_remaining	[expr $files_remaining - 2]
	set blast_file_list	[lremove $blast_file_list $file_A]
	set blast_file_list	[lremove $blast_file_list $file_B]
}