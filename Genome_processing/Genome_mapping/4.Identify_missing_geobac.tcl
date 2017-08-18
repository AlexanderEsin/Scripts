#!/usr/local/bin/tclsh
source ~/Dropbox/Scripts/General_utils.tcl
source ~/Dropbox/Scripts/Procs/gc_content.tcl
package require sqlite3

## Set paths ##
set direct /users/aesin/desktop
set working_directory $direct/Geo_analysis/Geo_ortholog_nucl

## Databases ##
sqlite3 all_prot_geo_db $direct/Clean_Proteomes/all_prot_geo_ncbi_db
sqlite3 genome_posit_db $working_directory/Geo_v_all_orth_location_db

cd $working_directory/Groups

set fasta_files [lsort -dictionary [glob *faa]]

set total_geobac 0
set total_present 0

set geobac_fasta_not_found {}
set geobac_table_not_found {}
set not_present_CDS_lengths {}
set present_CDS_lengths {}


## Set up the counter for total proteins ##
set tot_prots 0
foreach file $fasta_files {
	set prot_seqs [split_genes_fast [openfile $file]]
	set tot_prots [expr $tot_prots + [llength $prot_seqs]]
}

progress_init $tot_prots
set prots_processed 0


foreach file $fasta_files {
	set prot_seqs [split_genes_fast [openfile $file]]
	set group [string range $file 0 end-4]

	foreach prot_seq $prot_seqs {
		set prot_seq [string trim $prot_seq]
		set prot_id [string range $prot_seq 1 [string first " " $prot_seq]-1]
		
		set class [string range [all_prot_geo_db eval {SELECT class FROM t1 where id = $prot_id}] 1 end-1]
		
		if {$class eq "Geobac"} {
			set present [genome_posit_db eval {SELECT EXISTS(SELECT 1 FROM t1 WHERE prot_id = $prot_id LIMIT 1)}]
			set total_present [expr $total_present + $present]
			incr total_geobac

			if {$present == 0} {
				set header [string range $prot_seq 0 [string first \n $prot_seq]-1]
				set protein_name [string trim [string range $header [string first \) $header]+1 [string last \[ $header]-1]]

				set CDS [string trim [string range $prot_seq [string first \n $prot_seq] end]]
				regsub -all "\n" $CDS {} collapsed_CDS
				set CDS_length [string length $collapsed_CDS]

				lappend not_present_CDS_lengths $CDS_length
				lappend geobac_table_not_found "$group\t$prot_id\t$protein_name"
				lappend geobac_fasta_not_found $prot_seq
			} else {
				set CDS [string trim [string range $prot_seq [string first \n $prot_seq] end]]
				regsub -all "\n" $CDS {} collapsed_CDS
				set CDS_length [string length $collapsed_CDS]

				lappend present_CDS_lengths $CDS_length
			}
		}
		incr prots_processed
		progress_tick $prots_processed
	}
}

puts "Total Geobacilli queried: $total_geobac"
puts "Total with coordinates: $total_present"

## A table file with GROUP | PROT_ID | PROTEIN_NAME ##
set out [open $working_directory/Geobac_proteins_no_coords/Geobacillus_proteins_no_coordinates.tsv w]
puts $out [join $geobac_table_not_found \n]
close $out

## AA lengths of the proteins that cannot be mapped ##
set out [open $working_directory/Geobac_proteins_no_coords/Missing_AA_lengths.tsv w]
puts $out [join $not_present_CDS_lengths \n]
close $out

## AA lengths of the proteins that can be mapped ##
set out [open $working_directory/Geobac_proteins_no_coords/Present_AA_lengths.tsv w]
puts $out [join $present_CDS_lengths \n]
close $out

## A fasta file containing all the fasta sequences of the proteins that cannot be mapped ##

set out [open $working_directory/Geobac_proteins_no_coords/Geobac_not_mapped_fasta.faa w]
puts $out [join $geobac_fasta_not_found \n]
close $out

#puts "\n\n\n\n$total_geobac"
all_prot_geo_db close
genome_posit_db close


