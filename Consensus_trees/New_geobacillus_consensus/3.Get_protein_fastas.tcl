#!/usr/local/bin/tclsh

source "/Users/aesin/Documents/Scripts/General_utils.tcl"
package require sqlite3

## Open database to find the protein sequences based on protein IDs ##
sqlite3 db1 "/Users/aesin/Desktop/Clean_Proteomes/all_prot_geo_ncbi_db"

#############################################
## Read in necessary files and prep output ##
#############################################

## Read in list of vertical gene families which are monophyletic ##
set mono_vert_fams_file	"/Users/aesin/Desktop/Geo_analysis/New_geobacillus_consensus/Vertical_monophyletic_fams.tsv"
set mono_vert_fams_l	[split [string trim [openfile $mono_vert_fams_file]] \n]

puts stdout "There are [llength $mono_vert_fams_l] monophyletic vertical gene families"
# There are 738 monophyletic vertical gene families #

## Read in the table that has all the associated protein IDs (from 1.Get_vertical_fams_prot_ids.tcl) ##
set vert_prot_ids_file	"/Users/aesin/Desktop/Geo_analysis/New_geobacillus_consensus/Vertical_family_prot_ids.tsv"
set vert_prot_ids_l		[split [string trim [openfile $vert_prot_ids_file]] \n]

## Prepare output directory ##
set output_dir "/Users/aesin/Desktop/Geo_analysis/New_geobacillus_consensus/Gene_family_fastas"
file mkdir $output_dir

#########################################################
## Process the protein fastas for each relevant family ##
#########################################################

## For each vertical, monophyletic gene family get the list of protein IDs, get their sequence from the DB, and save to file ##
foreach mono_vert_fam $mono_vert_fams_l {
	
	set vert_prot_entry [lindex [lsearch -all -inline -glob $vert_prot_ids_l $mono_vert_fam\t*] 0]

	## Check that the appropriate entry has been found ##
	if {[string length $vert_prot_entry] == 0} {
		puts stdout "The prot id entry for gene family: $mono_vert_fam was not found. Error. Exiting..."
		exit
	}

	## Split the entry into individual protein IDs ##
	set vert_prot_ids [lrange [split $vert_prot_entry \t] 1 end]

	## For each protein ID - find corresponding fasta in table, add to per-family list ##
	set prot_id_per_family_l {}
	foreach protein_id $vert_prot_ids {
		set protein_fasta [lindex [db1 eval {select gene from t1 where id = $protein_id}] 0]

		# Check the the fasta has been found #
		if {[string length $protein_fasta] == 0} {
			puts stdout "The protein fasta entry for ID: $protein_id was not found. Error. Exiting..."
			exit
		}

		lappend prot_id_per_family_l $protein_fasta
	}

	## Check that all 22 protein fastas are present ##
	if {[llength $prot_id_per_family_l] != 22} {
		puts stdout "The number of protein fastas for family: $mono_vert_fam was not 22. Error. Exiting..."
		exit
	}

	## Process the fasta file and save it to necessary folder ##
	set prot_id_fasta_format [join $prot_id_per_family_l \n]
	set out [open $output_dir/Family_$mono_vert_fam\.faa w]
	puts $out $prot_id_fasta_format
	close $out

}

db1 close