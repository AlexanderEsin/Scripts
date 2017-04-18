#!/usr/local/bin/tclsh

## In this script we find gene families that show vertical descent into Geobacillus, and contain 1-to-1 orthologs in all Anoxybacillus and Geobacillus genomes ##
source ~/Documents/Scripts/General_utils.tcl
package require sqlite3


## Open database ##
sqlite3 db1 "/Users/aesin/Desktop/Geo_analysis/Geo_ortholog_nucl/Geo_v_all_orth_location_NEW_db"


## Get list of vertical gene families. Use 3,4,5,6 for the most conservative set ##
set vertical_set_file "/Users/aesin/Desktop/Mowgli/Consistent_HGT_Vertical/Consistent_Full/Scenarios_1_2/true_vertical_intersect_t3_t4_t5_t6.tsv"
set vertical_family_l [split [string trim [openfile $vertical_set_file]] \n]


## For each vertical gene family, get all Anoxy/Geobacillus binomials. Select where there are only 22, and all 22 are unique (1-to-1 orthologs) ##
set num_full_fams 0
set vertical_full_fam_l {}

set vertical_family_l1 [lrange $vertical_family_l 0 10]
foreach vertical_family $vertical_family_l {
	puts "Family: $vertical_family"

	# Query database. Lsearch faster than searching sqlite3 with glob binomial #
	set geobac_l_temp [db1 eval {select binomial from t1 where group_number = $vertical_family}]
	set geobac_l	[lsearch -all -inline -glob $geobac_l_temp "Geobacillus*"]
	set anoxybac_l	[lsearch -all -inline -glob $geobac_l_temp "Anoxybacillus*"]
	set anoxygeo_l	[concat $geobac_l $anoxybac_l]

	# Count the binomials #
	set num_anoxygeo [llength $anoxygeo_l]

	# Get unique #
	set anoxygeo_unique_l	[lsort -unique $anoxygeo_l]
	set num_anoxygeo_unique	[llength $anoxygeo_unique_l]

	puts "\tUnique Anoxyeo: $num_anoxygeo_unique"

	# Only 22, all 22 unique #
	if {$num_anoxygeo == 22 && $num_anoxygeo_unique == 22} {
		lappend vertical_full_fam_l $vertical_family
		incr num_full_fams
	}
}

puts $num_full_fams
# 738 gene Geobacillus families #


## For each gene family that has 22 vertically descended 1-to-1 orthologs, get the protein IDs. Ensure all 22 IDs are found ##
set output_tbl_l {}
foreach vertical_full_fam $vertical_full_fam_l {
	set geo_prot_id_l	[db1 eval {select prot_id from t1 where group_number = $vertical_full_fam and binomial glob "Geobacillus*"}]
	set anoxy_prot_id_l	[db1 eval {select prot_id from t1 where group_number = $vertical_full_fam and binomial glob "Anoxybacillus*"}]
	
	set anoxygeo_prot_id_l [concat $geo_prot_id_l $anoxy_prot_id_l]
	set num_prot_id [llength $anoxygeo_prot_id_l]

	puts "Family: $vertical_full_fam\tProtein IDs: $num_prot_id"
	if {$num_prot_id != 22} {
		puts "Not 22 proteins. Error"
		exit
	}

	set fam_prot_entry [join [concat $vertical_full_fam $anoxygeo_prot_id_l] \t]
	lappend output_tbl_l $fam_prot_entry

}

## Write out the results in format: Gene family \t prot_id_1 \t prot_id_2 \t etc ... ##
set out [open "/Users/aesin/Desktop/Geo_analysis/New_geobacillus_consensus/Vertical_family_prot_ids.tsv" w]
puts $out [join $output_tbl_l \n]
close $out


# Close DB #
db1 close