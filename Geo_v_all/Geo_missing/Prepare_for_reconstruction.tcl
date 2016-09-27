#!/usr/local/bin/tclsh

### This script takes all the Groups that had lost some Geobacillus sequences during the context script that picked a better level at which to reconsutruct the family (see Missing_geobac_families.tsv) and prepares them for reconstruction at a level where all Geobacillus sequences will be included. Total number of Geobacillus sequences that were "lost" = 1,304, belonging to 283 groups ###

source ~/Dropbox/Scripts/General_utils.tcl
package require sqlite3

###########################################################################
## Set variables ##
set direct /users/aesin/desktop/Geo_analysis/Geo_v_all/2.0
set missing_geo_path $direct/Missing_geobac_reconstruction
set final_group_fastas_path $missing_geo_path/Final_family_groups_all_geo

###########################################################################
## Get the list of missing families ##
set missing_geo_data [lrange [split [string trim [openfile $missing_geo_path/Missing_geobac_families.tsv]] \n] 1 end]
set missing_geo_l {}
foreach group $missing_geo_data {lappend missing_geo_l [lindex [split $group \t] 0]}

###########################################################################
## For each group that had missing Geobacillus sequences, get the fasta at -10 ##
file mkdir $final_group_fastas_path
set total_prot_num 0

foreach group $missing_geo_l {
	file copy -force $final_group_fastas_path/$group\.faa $missing_geo_path/Fastas
	set total_prot_num [expr $total_prot_num + [llength [split_genes [openfile $final_group_fastas_path/$group\.faa]]]]
}

puts "Total number of proteins to be processed from [llength $missing_geo_l] groups: $total_prot_num"
## Check how many proteins we will have to process in this set ##