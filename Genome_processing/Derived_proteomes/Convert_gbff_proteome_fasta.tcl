#!/usr/local/bin/tclsh

## In this script we create a set of proteomes based on the genbank files in Clean_Genomes. These can then be used as a reference to blast a sequence against to find an updated protein id ##

source ~/Dropbox/Scripts/General_utils.tcl

###########################################################################

set direct /users/aesin/desktop/Clean_Genomes
set gbff_file_dir Prok_gbff_files
set output_dir Gbff_derived_proteomes

file mkdir $direct/$output_dir

cd $direct/$gbff_file_dir
set gbff_files [glob *gbff]
set counter 1

foreach gbff_file $gbff_files  {
	set output_name "[string range $gbff_file 0 end-5].faa"
	if {[file exists $direct/$output_dir/$output_name] == 1} {
		puts "Already converted... skipping"
		continue
	} else {
		catch {exec ADE_genbank_to_fasta.py -i $gbff_file -o $direct/$output_dir/$output_name -d tab -q "protein_id,product,locus_tag,accessions,location"}
		puts "Converted gbff --> proteome fasta --- $counter / [llength $gbff_files]"
	}
	incr counter
}

puts "Files to be converted: [llength $gbff_files]"
puts "Total files after conversion: [llength [glob $direct/$output_dir/*faa]]"

