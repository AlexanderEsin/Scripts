#!/usr/local/bin/tclsh

source ~/Documents/Scripts/General_utils.tcl


set direct /users/aesin/Desktop/Geo_again/AGeo_inparalogs

## Input and DB copied from the cross-Anogeo blast
set input_path	$direct/Input
set db_path		$direct/DB
set out_path	$direct/Blast_out

## Lists of query proteomes and corresponding DBs
cd $input_path
set proteome_l [glob *.fasta]
cd $db_path
set db_l [glob *.pin]

## Output path
file mkdir $out_path

## Foreach proteome run blast against itself
## Run the blast at 1e-50 - don't expect reliable in-paralogs at lower thresholds
## Using desktop blast: v2.6.0+
foreach query $proteome_l {
	set query_base	[string range $query 0 end-15]
	set db_base		$query_base

	catch {exec blastp -query $input_path/$query -db $db_path/$db_base -out $out_path/$query_base.tsv -evalue 1e-50 -outfmt 6 -max_hsps 1 -num_threads 2}
}

## Take the Blast output files and trim away blast hits 
## against the query protein (self hits)
file mkdir $direct/Self_trimmed
file mkdir $direct/Poss_paralogs

set paralog_stats	{"Proteome\tNum.Possible.Paralogs"}

cd $out_path
set blast_out_files	[glob *tsv]
foreach blast_output $blast_out_files {

	## Get accession_assembly and write a ticker
	set acc_ass	[string range $blast_output 0 end-4]
	puts "Finding paralogs in $blast_output"

	## Trim away self hits
	catch {exec /Users/aesin/Documents/Scripts/Geo_again/AGeo_inparalogs/1.Trim_selfhits.R $out_path/$blast_output $direct/Self_trimmed}

	## File name remains the same in the new directory
	## Do a trivial RBBH - keeping bitscores
	set hits	[split [string trim [openfile $direct/Self_trimmed/$blast_output]] \n]

	set paralog_counter 0
	set paralog_list	{}

	foreach blast_hit $hits {
		set cols	[split $blast_hit \t]
		set gene_1	[lindex $cols 0]
		set gene_2	[lindex $cols 1]

		## Make sure the reverse has not already been done!
		if {[lsearch $paralog_list "$gene_2\t$gene_1\t*"] != -1} {
			continue
		}

		set RBBH_hit	[lsearch -all -inline $hits "$gene_2\t$gene_1*"]

		if {[llength $RBBH_hit] == 1} {
			set hit		[lindex $RBBH_hit 0]
			set bits_1	[lindex $cols 11]
			set bits_2	[lindex [split $hit \t] 11]

			set eval_1	[lindex $cols 10]
			set eval_2	[lindex [split $hit \t] 10]

			set mean_bits	[expr ($bits_1 + $bits_2) / 2]
			set mean_eval	[expr ($eval_1 + $eval_2) / 2]
		} else {
			continue
		}

		## Only writing out 1 line for both genes, so when searching for paralogs
		## we must check both columns
		lappend paralog_list "$gene_1\t$gene_2\t$mean_bits\t$mean_eval"
		incr paralog_counter
	}

	set out [open $direct/Poss_paralogs/$blast_output w]
	puts $out [join $paralog_list \n]
	close $out

	## Write the number of possible paralogs into a table
	lappend paralog_stats "$acc_ass\t$paralog_counter"
}


## Write out the number of potential paralogs
## detected per genome
set out [open $direct/Possible_paralog_numbers.tsv w]
puts $out [join $paralog_stats \n]
close $out


## // ##

## Write out a global Anogeo paralog table for use downstream for comparisons
cd $direct/Poss_paralogs
set paralog_files	[glob *tsv]

## For each per-genome paralog file, open and
## combine into a global file
set anogeo_paralog_tbl	{}
foreach paralog_file $paralog_files {
	set paralog_data [string trim [openfile $paralog_file]]
	lappend anogeo_paralog_tbl $paralog_data
}

## Write out a global pairwise table of possible paralogs
set out [open $direct/Anogeo_possParalog_tbl.tsv w]
puts $out [join $anogeo_paralog_tbl \n]
close $out









