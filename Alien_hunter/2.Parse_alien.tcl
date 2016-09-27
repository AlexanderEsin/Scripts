#!/usr/local/bin/tclsh

###########################################################################
## Procs ##

source ~/Dropbox/Scripts/General_utils.tcl
package require sqlite3

###########################################################################

set direct /users/aesin/desktop/Parametric_HGT

cd /users/aesin/Desktop/Geo_analysis/Geo_ortholog_nucl
sqlite3 ortholog_position_db Geo_v_all_orth_location_db

###########################################################################

set input_dir $direct/Alient_hunter_output
set output_dir $direct/Predicted_genes
file mkdir $output_dir

# puts "Number of genomes: [llength $geo_genomes]"

cd $input_dir
set output_dirs [glob -type d *]
set dir [lindex $output_dirs 0]

#foreach dir $output_dirs {
	set ncbi_id $dir

	cd $input_dir/$dir
	set output_data [string trim [openfile $ncbi_id\_prediction.txt]]

	set hgt_hits [regexp -all -line -inline {.+(?:misc_feature).+} $output_data]



	#puts $ncbi_id
	file mkdir $output_dir/$ncbi_id

	#puts $output_dir/$ncbi_id/$ncbi_id\_prediction.txt

	$alien_hunter_path/alien_hunter $genome $output_dir/$ncbi_id/$ncbi_id\_prediction.txt
#}

ortholog_position_db close