#!/usr/local/bin/tclsh
source ~/Dropbox/Scripts/General_utils.tcl
source ~/Dropbox/Scripts/Procs/gc_content.tcl

## Set variables ##
set folder Prok_genomic
set direct ~/desktop/Clean_genomes

## Make output folder ##
file mkdir $direct/$folder\_fna_reform

## Get list of fna files ##
cd $direct/$folder\_fna
set unformatted_files [lsort -dictionary [glob *fna]]

## Set empty variables and counter ##
set fna_header_db {}
set element_ids {}
set i 1

## For each file, reformat the fna (see ~/Dropbox/Scripts/Procs/gc_content.tcl) and write out the reformatted file ##
foreach fna_file $unformatted_files {

	reform_fna_fasta $fna_file

	lappend fna_header_db $header_db_entry
	lappend element_ids $element_id_list

	set out [open $direct/$folder\_fna_reform/$fna_file w]
	puts $out $no_gaps_fasta
	close $out

	puts "Reformatted .. $i / [llength $unformatted_files]. Total entries in the genome element DB: [llength $fna_header_db]"
	incr i
}

set element_ids [join $element_ids \n]
set element_id_split [split $element_ids \n]

## Check whether any element ids (e.g. chromosome accessions) have duplicates - typically this should NOT happen ##
set duplicated_element_ids [dups $element_id_split]
if {[llength $duplicated_element_ids] > 0} {
	puts "The following element ids have duplicates (check the $folder\_genome_element_DB.tsv file):\n[join $duplicated_element_ids \n]"
} else {
	puts "No element ids are duplicated."
}

## Comprehensive database of all the elements parsed ##
set out [open $direct/$folder\_genome_element_DB.tsv w]
puts $out [join $fna_header_db \n]
close $out