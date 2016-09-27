#!/usr/local/bin/tclsh

### Get mean lengths of Geobacillus proteins ###

source ~/Dropbox/Scripts/General_utils.tcl
package require math::statistics

set direct /users/aesin/desktop/Geo_analysis/Geo_v_all/2.0/Geo_AG_only_groups
file mkdir $direct/Stats

cd $direct/Combined_map_fastas
set group_fastas [lsort -dictionary [glob *faa]]
set i 1

set output_table [list "Gene_number\tMean"]

foreach fasta $group_fastas {
	set group_number [string range $fasta 0 end-4]
	set genes [split_genes_fast [openfile $fasta]]
	set CDS_lengths {}

	foreach gene $genes {
		if {[regexp "Geobacillus" $gene] == 1} {
			regsub -all "\n" [string range $gene [string first \n $gene] end] {} CDS
			set CDS_len [string length $CDS]
			lappend CDS_lengths $CDS_len
		}
	}

	set mean_geo_len [::math::statistics::mean $CDS_lengths]
	# set len_stdev [::math::statistics::stdev $CDS_lengths]
	lappend output_table "[llength $CDS_lengths]\t$mean_geo_len"
	puts "$i / [llength $group_fastas]"
	# if {$i != $group_number} {
	# 	puts "ERROR $group_number"
	# 	exit
	# }
	incr i
}

set out [open $direct/Stats/Combined_map_protein_size_mean.tsv w]
puts $out [join $output_table \n]
close $out