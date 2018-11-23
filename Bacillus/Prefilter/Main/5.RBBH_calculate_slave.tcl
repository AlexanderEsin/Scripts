#!/usr/bin/tclsh

source /home/ade110/Scripts/General_utils.tcl

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
set direct		/home/ade110/Work/Bacillus/Prefilter/Bacillus_withOut
set blast_split	$direct/BLAST/OUT_all_rearrange

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
		exit 1
	}

	# Check that it's not a self blast - if it is, skip
	if {$file_A eq $file_B} {
		set files_remaining	[expr $files_remaining - 1]
		set blast_file_list	[lremove $blast_file_list $file_A]
		puts "Self blast -- skipping..."
		continue
	}

	set rbbh_output			{"GeneA\tGeneB\tEval"}
	set paralog_comparison	{"GeneA\tGeneB\tBits"}
	set rbbh_counter 0
	set below_cutoff 0
	set miss_counter 0
	set hits_processed 0

	set hits_A	[split [string trim [openfile $file_A]] \n]
	set hits_B	[split [string trim [openfile $file_B]] \n]

	### For some reason the max_hsps_per_subject option did not work
	### So, we need to select only the best alignment for each blast
	### (it seems max_subject did work - so we shouldn't have multiple subjects)

	## First, select all the unique query genes in A
	set unique_genes_A	{}
	foreach blast_hit $hits_A {
		set cols_A		[split $blast_hit \t]
		set gene_1		[lindex $cols_A 0]
		lappend unique_genes_A $gene_1
	}
	set unique_genes_A	[lsort -unique $unique_genes_A]

	## For each unique gene A, select the best blast alignment with B
	## This will be the first occurence of this genes, as output is ordered
	## by evalue
	set best_hits_A		{}
	foreach unique_gene_A $unique_genes_A {
		set all_hits_A	[lsearch -all -inline $hits_A "$unique_gene_A\t*"]
		set best_hit_A	[lindex $all_hits_A 0]
		lappend best_hits_A	$best_hit_A
	}


	## Now, for each best alignment in A, we find the best alignment in B
	## and ensure that the gene IDs are reciprocal
	foreach blast_hit $best_hits_A {
		incr hits_processed

		set cols_A	[split $blast_hit \t]
		set gene_1	[lindex $cols_A 0]
		set gene_2	[lindex $cols_A 1]

		## For each best blast of a gene in A against B, check wether 
		## the subject gene (in B) has any blast hits against A
		set all_hits_B	[lsearch -all -inline $hits_B "$gene_2\t*"]
		if {[llength $all_hits_B] == 0} {
			incr miss_counter
			continue
		} else {
			set best_hit_B	[lindex $all_hits_B 0]
		}

		## Is the best hit of gene 2 in B gene 1 in A? (RBBH)
		## If not, this is a miss
		set best_hit_B_subject	[lindex [split $best_hit_B \t] 1]
		if {$best_hit_B_subject eq $gene_1} {
			set RBBH_hit	$best_hit_B
		} else {
			incr miss_counter
			continue
		}


		## If this is a true RBBH, proceed to extract evalues
		set eval_A	[lindex $cols_A 10]
		set eval_B	[lindex [split $RBBH_hit \t] 10]

		if {$eval_cutoff == 1e-10} {

			## Only write out the paralog comparison at eval = 1e-10
			set bits_A	[lindex $cols_A 11]
			set bits_B	[lindex [split $RBBH_hit \t] 11]
			set av_bits	[expr ($bits_A + $bits_B) / 2]

			## A condensed format for in-paralog comparison using bitscore
			lappend paralog_comparison	"$gene_1\t$gene_2\t$av_bits"

			## Append the full RBBH output
			lappend rbbh_output "$gene_1\t$gene_2\t$eval_A"
			lappend rbbh_output "$gene_2\t$gene_1\t$eval_B"
			incr rbbh_counter

		## If the evalue cutoff is not at default minimum for blast (1e-10)
		## check that it's below the treshold both ways
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
	}

	## Prepare the output file name for the RBBH hits
	set RBBH_out_name	[string range $file_A 0 end-4]

	## Write out the RBBH hits
	set out		[open $RBBH_out/$RBBH_out_name\.txt w]
	puts $out	[join $rbbh_output \n]
	close $out

	## Write out the condensed RBBH with bitscores for comparison
	## to in-paralogs. Only at 1e-10 (all hits)
	if {$eval_cutoff == 1e-10} {
		set out		[open $para_out/$RBBH_out_name\.txt w]
		puts $out	[join $paralog_comparison \n]
		close $out
	}

	## Counter for loop
	set files_remaining	[expr $files_remaining - 2]
	set blast_file_list	[lremove $blast_file_list $file_A]
	set blast_file_list	[lremove $blast_file_list $file_B]
}