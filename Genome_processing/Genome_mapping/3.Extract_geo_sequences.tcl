#!/usr/local/bin/tclsh
source ~/Dropbox/Scripts/General_utils.tcl
package require sqlite3

## Set paths ##
set direct /users/aesin/desktop
set working_directory $direct/Geo_analysis/Geo_ortholog_nucl
set output_map_directory $working_directory/Groups_all_map_geo_only
set output_unmap_directory $working_directory/Groups_not_all_map_geo_only

file mkdir $output_map_directory
file mkdir $output_unmap_directory

## Databases ##
sqlite3 genome_posit_db $working_directory/Geo_v_all_orth_location_NEW_db


###########################################################################
## From each of the groups that were still eligible after removing any Geobacillus sequences that did not map to the new annotation, extract all the Geobacillus proteins sequences remaining into separate files ##
cd $working_directory/Groups_all_geo_map
set groups [lsort -dictionary [glob *faa]]

set total_geo_processed 0

foreach group $groups {
	set geo_prot_out_l {}
	set prots [split_genes [openfile $group]]
	foreach prot $prots {
		set prot_id [string trim [string range $prot 1 [string first " " $prot]]]
		set class [genome_posit_db eval {SELECT class FROM t1 WHERE prot_id = $prot_id}]

		if {$class eq "Geobac"} {
			lappend geo_prot_out_l $prot
			incr total_geo_processed
		}
	}
	set out [open $output_map_directory/$group w]
	puts $out [join $geo_prot_out_l {}]
	close $out
	
}

puts "Total geobacillus extracted: $total_geo_processed"
genome_posit_db close

###########################################################################
## Extract all the geobacillus from the families that did not have all the Geobacillus map to genome coordinates (for ontology testing) ##

sqlite3 all_prot_geo_db $direct/Clean_Proteomes/all_prot_geo_ncbi_db

cd $working_directory/Groups_not_all_geo_map
set groups [lsort -dictionary [glob *faa]]

set total_geo_processed 0

foreach group $groups {
	set geo_prot_out_l {}
	set group_data [openfile $group]
	set prots [split_genes $group_data]

	foreach prot $prots {
		set prot_id [string trim [string range $prot 1 [string first " " $prot]]]
		set class ""
		set class [all_prot_geo_db eval {SELECT class FROM t1 WHERE id = $prot_id}]
		## In all_prot_geo_ncbi_db the class is surrounded by curly brackets ##
		set class [string range $class 1 end-1]
		if {$class eq ""} {
			puts stderr "ERROR: $group\n$prot"
			exit
		}

		if {$class eq "Geobac"} {
			lappend geo_prot_out_l $prot
			incr total_geo_processed
		}
	}
	set out [open $output_unmap_directory/$group w]
	puts $out [join $geo_prot_out_l {}]
	close $out
	
}

puts "Total geobacillus extracted: $total_geo_processed"
all_prot_geo_db close