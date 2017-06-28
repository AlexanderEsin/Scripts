#!/usr/local/bin/tclsh

source ~/Documents/Scripts/General_utils.tcl
package require sqlite3

# Open sqlite database
sqlite3 db1 /Users/aesin/Desktop/Clean_Proteomes/all_prot_geo_ncbi_db

## 1. Get list of inside taxa ##
set fasta_file_path /Users/aesin/Desktop/Clean_Proteomes/Inside_group

cd $fasta_file_path
set in_group_fastas [glob -nocomplain *faa]

# Unzip if necessary - if unzipped, then rezip afterwards
if {[llength $in_group_fastas] == 0} {
	set zipped_fastas [glob -nocomplain *faa.gz]
	
	if {[llength $zipped_fastas] == 0} {
		puts "Cannot find fasta or zipped fasta files in $fasta_file_path"
		exit 1
	} else {
		set rezip TRUE
		foreach zip_file $zipped_fastas {exec gunzip $zip_file}

		set in_group_fastas [glob -nocomplain *faa]

		if {[llength $in_group_fastas] > 0} {
			puts "Fasta files unzipped"
		} else {
			puts "Unzipping failed"
		}
	}
} else {
	set rezip FALSE
}


## 2. Get binomial name for each file ##
set binomial_l {}
foreach fasta_file $in_group_fastas {
	set binomial [db1 eval {SELECT DISTINCT binomial FROM t1 WHERE file_name = $fasta_file}]
	lappend binomial_l $binomial
}
puts "Total binomials belonging to `Inside Group`: [llength $binomial_l]"
# Close DB
db1 close


## 3. (OPTIONAL) rezip input fasta files ##
if {$rezip == TRUE} {
	foreach fasta_file $in_group_fastas {exec gzip $fasta_file}
}


## 4. Confirm all IG group binomials are present in mowgli species tree ##
set tree_data [openfile /Users/aesin/Desktop/Mowgli/Species_tree/Ultrametric_species_tree_tipfix.txt]

foreach binomial $binomial_l {
	if {[regexp $binomial $tree_data] == 0} {
		puts "Cannot find $binomial in the species tree"
		exit 1
	}
}

file mkdir /Users/aesin/Desktop/Mowgli/Inside_group
set out [open /Users/aesin/Desktop/Mowgli/Inside_group/Inside_group_tips.txt w]
puts $out [join $binomial_l \n]
close $out