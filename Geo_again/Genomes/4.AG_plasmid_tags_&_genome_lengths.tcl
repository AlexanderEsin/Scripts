#!/usr/local/bin/tclsh

## In this script we determine which locus IDs correspond to
## plasmid-bourned proteins. We only need to do this for the 
## AG genomes

source ~/Documents/Scripts/General_utils.tcl

## Directories
set direct		/users/aesin/Desktop/Geo_again/Genomes
set gbff_dir	$direct/Genome_gbffs
set plasmid_out	$direct/Plasmid_per_genome
file mkdir		$plasmid_out


## List of AG acc_ass
set AG_acc_asses	[split [string trim [openfile $direct/Genome_lists/AG_acc_Ass_names.txt]] \n]

## Get the input gbff files
set input_zipped		[glob -nocomplain $gbff_dir/*gz]
if {[llength $input_zipped] != 0} {
	puts stdout "Unzipping gbff files ..."
	catch {exec gunzip -r -f $gbff_dir}
	puts stdout "Unzipping gbff files ... DONE"
}
set input_gbffs		[glob $gbff_dir/*gbff]
set num_input_gbffs	[llength $input_gbffs]

## Pattern to split the gbk file into partitions 
## (e.g. either contigs or chromosome/plasmid)
set partition_pat {\n+(//)+\n+}

## Set output lists
set genome_len_tbl		{"Genome\tGenome_length_ex_plasmid"}
set AG_global_plasmids	{}
set counter				0

foreach gbff_file $input_gbffs {

	## Get the Acc_Ass from the file name
	## The _genomic.gbff suffix is 13 character
	set file_name	[file tail $gbff_file]
	set acc_ass		[string range $file_name 0 end-13]

	puts stdout "Processing $acc_ass - $counter // $num_input_gbffs"

	## See whether this genome is an AG genome
	if {[lsearch $AG_acc_asses $acc_ass] != -1} {
		set ag_genome	1
	} else {
		set ag_genome	0
	}

	## Get the gbff genome data and divide by partition (chromosome)
	## We can't split directly by the partition pattern, 
	## so we do an intermediate regsub
	set genome_data [string trim [openfile $gbff_file]]
	if {[regexp {£} $genome_data] == 1} {
		error "Gbff file should not have the \"£\" character natively"
	}
	regsub -all $partition_pat $genome_data "\n£\n" genome_data
	set partitions [split $genome_data {£}]


	## Counter and holders
	set total_length 0
	set total_clean_tags {}
	set plasmid_stats_file {}

	## Process each partition in turn
	foreach partition $partitions {

		## Check we have a valid partition
		if {[string length $partition] == 0} {
			error "Not expecting 0-length partitions"
		}
		set partition [string trim $partition]

		## For each partition extract the "LOCUS" annotation
		## This contains the chromosome accession and length;
		## make it tab-separated
		set LOCUS_line [string range $partition 0 [string first \n $partition]]
		regsub -all {\s+} $LOCUS_line "\t" LOCUS_line


		## Trim away the LOCUS line and get the next line
		## This is the "DEFINITION" annotation and contains the name
		## of the DNA element (including whether it's a plasmid)
		set partition_trim	[string trim [string range $partition [string first \n $partition] end]]
		set DEFINITION_line	[string range $partition_trim 0 [string first \n $partition_trim]]

		set plasmid_presence	[regexp {plasmid} $DEFINITION_line]

		## If this partition represents a plasmid, prepare the output file
		## Per-genome this is a list of partitions and associated locus tags
		if {$plasmid_presence == 1} {

			## Partition information
			lappend plasmid_stats_file $LOCUS_line
			lappend plasmid_stats_file $DEFINITION_line

			## All locus tags in this partition
			set locus_tags [regexp -all -inline -line {(?:/locus_tag=)+.+} $partition]
			set clean_partition_tags {}
			
			## String process the locus tags to get the string only
			foreach tag $locus_tags {
				set tag_trim [string range $tag [string first "\"" $tag]+1 [string last "\"" $tag]-1]
				lappend clean_partition_tags $tag_trim
			}

			## Remove any duplicates - though we don't expect any
			set clean_partition_tags	[lsort -unique $clean_partition_tags]
			set plasmid_stats_file		[concat $plasmid_stats_file	$clean_partition_tags]
			set total_clean_tags		[concat $total_clean_tags $clean_partition_tags]
		} else {
			set total_length [expr $total_length + [lindex [split $LOCUS_line \t] 2]]
		}
	}

	## Add the total length of the genome to the global table
	## Also add the plasmid locus tags to a global list
	lappend genome_len_tbl "$acc_ass\t$total_length"
	set AG_global_plasmids	[concat $AG_global_plasmids $total_clean_tags]

	## Write out the plasmid tags
	if {[llength $total_clean_tags] > 0} {
		set out [open $plasmid_out/$acc_ass\_plasmid_tags.tsv w]
		puts $out [join $total_clean_tags \n]
		close $out

		set out [open $plasmid_out/$acc_ass\_plasmid_stats.tsv w]
		puts $out [string trim [join $plasmid_stats_file \n]]
		close $out
	}
	incr counter
}

## Write out the genome lengths 
set out [open $direct/Genome_lengths.tsv w]
puts $out [join $genome_len_tbl \n]
close $out

## Write out the global plasmid tag table
set out [open $direct/All_plasmid_tags.tsv w]
puts $out [join $AG_global_plasmids \n]
close $out


