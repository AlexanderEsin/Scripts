#!/usr/local/bin/tclsh

## In this script we label each protein ID in each proteome with
## a unique suffix corresponding to the genome taxid and start position

## In addition, we remove any protein entries with a "missing protein ID"

source ~/Documents/Scripts/General_utils.tcl

set direct		/users/aesin/Desktop/Geo_again

set all_in_dir	$direct/BLAST/IN_all
set all_db_dir	$direct/BLAST/DB_all

set geo_in_dir	$direct/BLAST/IN_geo
set geo_db_dir	$direct/BLAST/DB_geo

## List of AG acc_ass
set AG_acc_asses	[split [string trim [openfile $direct/Genomes/Genome_lists/AG_acc_Ass_names.txt]] \n]


## For AG proteome, copy the proteome (input) and DB into respective directories
foreach AG $AG_acc_asses {
	set proteome	$all_in_dir/$AG\_proteome.fasta
	set db			[glob $all_db_dir/$AG*]
	
	if {[file exists $proteome] == -1 || [llength $db] < 3} {
		error "Can't find AG proteome: $AG"
	}

	file copy -force $proteome $geo_in_dir/
	foreach db_file $db {
		file copy -force $db_file $geo_db_dir/
	}
}