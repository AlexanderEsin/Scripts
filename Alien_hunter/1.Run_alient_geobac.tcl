#!/usr/local/bin/tclsh

###########################################################################
## Procs ##

source ~/Dropbox/Scripts/General_utils.tcl

###########################################################################

set direct /users/aesin/desktop
set alien_hunter_path $direct/Software/alien_hunter

set genomes_path $direct/Geo_analysis/Geo_omes/Geobac_genomes_raw

set output_dir $direct/Parametric_HGT/Alient_hunter_output
file mkdir $output_dir

cd $genomes_path
set geo_genomes [glob *fna]

puts "Number of genomes: [llength $geo_genomes]"

foreach genome $geo_genomes {
	set ncbi_id [string range $genome 0 [string last "_" $genome]-1]
	#puts $ncbi_id
	file mkdir $output_dir/$ncbi_id

	#puts $output_dir/$ncbi_id/$ncbi_id\_prediction.txt

	$alien_hunter_path/alien_hunter $genome $output_dir/$ncbi_id/$ncbi_id\_prediction.txt
}