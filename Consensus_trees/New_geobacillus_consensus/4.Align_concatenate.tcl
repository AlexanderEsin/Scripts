#!/usr/local/bin/tclsh

source "/Users/aesin/Documents/Scripts/General_utils.tcl"
package require sqlite3
package require math::statistics

###############################
## Paths and I/O directories ##
###############################

## Open database to find the protein sequences based on protein IDs ##
sqlite3 db1 "/Users/aesin/Desktop/Clean_Proteomes/all_prot_geo_ncbi_db"

## Input directory - contains fasta files per gene family ##
set raw_fasta_dir "/Users/aesin/Desktop/Geo_analysis/New_geobacillus_consensus/Gene_family_fastas"

## Output directory for the MSA files ##
set muscle_dir "/Users/aesin/Desktop/Geo_analysis/New_geobacillus_consensus/Gene_family_MSA"
file mkdir $muscle_dir

## Output directory for the per-species concatenated alignments ##
set per_species_align_dir "/Users/aesin/Desktop/Geo_analysis/New_geobacillus_consensus/Per_species_concat"
file mkdir $per_species_align_dir

####################################
## Run MUSCLE on each gene family ##
####################################

## Get input file list ##
set input_files [lsort -dictionary [glob -type f $raw_fasta_dir/*]]

#set input_files [lrange $input_files 0 100]; ## temp
set i 0
set bad_length 0

## Run muscle for each fasta file and output to MSA directory ##
foreach fasta_file $input_files {
	set exclude FALSE
	## Extract just file name without exension. E.g. Family_100 ##
	set raw_file_name [string range $fasta_file [string last \/ $fasta_file]+1 end-4]

	set proteins [split_genes [openfile $fasta_file]]
	set prot_length_l {}
	foreach protein $proteins {
		regsub -all "\n" [string trim [string range $protein [string first \n $protein] end]] {} seq_no_gaps
		set seq_len [string length $seq_no_gaps]
		lappend prot_length_l $seq_len
	}
	set len_mean [::math::statistics::mean $prot_length_l]
	#puts stdout "Mean for family: $raw_file_name is $len_mean"
	set quarter_mean [expr $len_mean / 4]
	set outlier_max [expr $len_mean + $quarter_mean]
	set outlier_min [expr $len_mean - $quarter_mean]

	foreach length $prot_length_l {
		if {$length < $outlier_min || $length > $outlier_max} {
			puts "In family $raw_file_name there is a protein with an outlier length"
			set exclude TRUE
		}
	}

	## Run muscle ##
	incr i
	if {$exclude == FALSE} {
		catch {exec muscle -in $fasta_file -out $muscle_dir/$raw_file_name\.fas}
		puts stdout "--> MUSCLE $i / [llength $input_files]"
	} else {
		puts stdout "--> SKIP: Protein length in family inconsistent"
		incr bad_length
	}
}

puts "\n\n$bad_length / $i gene families that did not conform to the consistent length requirements"

####################################
## Stitch the alignments together ##
####################################

set msa_files [lsort -dictionary [glob -type f $muscle_dir/*]]

## Set up output files based on binomial names. Use the first MSA to derive binomial names to be used for the rest. ##
set first_msa_file [lindex $msa_files 0]
set msa_proteins [split_genes [openfile $first_msa_file]]

set open_file_channels {}
foreach protein $msa_proteins {
	set protein_id [string range $protein 1 [string first " " $protein]-1]
	set binomial [db1 eval {select binomial from t1 where id = $protein_id}]
	## Open output files - one for each binomial, based on the binomial name ##
	set $binomial\_out [open $per_species_align_dir/$binomial\.fas w+]
	lappend open_file_channels $binomial\_out
}

## For each alignment, extract the sequences and place them in the relevant per-binomial file ##
set total_alignment_length 0
set family_position_l {}
foreach msa_file $msa_files {
	set msa_proteins [split_genes [openfile $msa_file]]

	## For each protein get the binomial (to match to file), and also the flat CDS - i.e. no newline breaks ##
	foreach protein $msa_proteins {
		set protein_id [string range $protein 1 [string first " " $protein]-1]
		set binomial [db1 eval {select binomial from t1 where id = $protein_id}]
		
		regsub -all \n [string trim [string range $protein [string first \n $protein] end]] {} flat_CDS
		
		## Add the flat CDS to the correct per-binomial concatenated file ##
		puts [expr $$binomial\_out] $flat_CDS
	}
	
	set family_num [string range $msa_file 0 end-4]
	set fam_pos_entry "$family_num\t$total_alignment_length"
	lappend family_position_l $fam_pos_entry

	## All flat CDS should be same length ##
	set total_alignment_length [expr $total_alignment_length + [string length $flat_CDS]]
}

## Close all the open per-binomial alignment files##
foreach open_file $open_file_channels {
	close [expr $$open_file]
}

set out [open "/Users/aesin/Desktop/Geo_analysis/New_geobacillus_consensus/Family_position_in_concat.tsv" w]
puts $out [join $family_position_l \n]
close $out

puts "Total alignment length written out all of the individual files: $total_alignment_length"

########################################################################
## Combine individual concatenated alignments into a master alignment ##
########################################################################

set per_species_concat_files [glob -type f $per_species_align_dir/*]
set master_alignment_l {}

## Open and remove all newlines generated by the writing out process in the per-binomial concat files ##
foreach per_species_concat_file $per_species_concat_files {
	set newlines_removed [regsub -all "\n" [string trim [openfile $per_species_concat_file]] {} no_newlines_data]
	set out [open $per_species_concat_file w]
	puts $out $no_newlines_data
	close $out
	puts "For file $per_species_concat_file\: $newlines_removed newlines removed"
}

## For each per-species file, extract sequence - linebreak it into 80 characters per line, add header (binomial), and add to master alignment file ##
set key_l {}
set k 1
foreach per_species_concat_file $per_species_concat_files {
	set binomial_name [string range $per_species_concat_file [string last \/ $fasta_file]+1 end-4]

	set concat_sequence [string trim [openfile $per_species_concat_file]]
	set linebroken_concat_seq [linebreak $concat_sequence]

	set fasta_format_seq ">k$k\n$linebroken_concat_seq"
	lappend master_alignment_l $fasta_format_seq

	lappend key_l "k$k\t$binomial_name"
	incr k
}

set out [open "/Users/aesin/Desktop/Geo_analysis/New_geobacillus_consensus/Master_concat_alignment.fas" w]
puts $out [join $master_alignment_l \n]
close $out

set out [open "/Users/aesin/Desktop/Geo_analysis/New_geobacillus_consensus/Master_alignment_key.txt" w]
puts $out [join $key_l \n]
close $out

exec /bin/bash -c "seqret '/Users/aesin/Desktop/Geo_analysis/New_geobacillus_consensus/Master_concat_alignment.fas' '/Users/aesin/Desktop/Geo_analysis/New_geobacillus_consensus/Master_concat_alignment.phy' -osformat phylip 2> /dev/null"

db1 close