#!/usr/local/bin/tclsh
source ~/Documents/Scripts/General_utils.tcl
package require sqlite3

set master_dir		/Users/aesin/Desktop/Geo_again
set output_dir		$master_dir/Functional_annotation/AG_only_groups
set famGrp_dir		$master_dir/Family_mapping/Reduced_groups_fasta
set protDB_path		$master_dir/All_prot_db

file mkdir $output_dir

# Open database
sqlite3 all_prot_db $protDB_path

# List of reduced fasta groups (input to tree building)
set full_groups_files	[glob $famGrp_dir/*fasta]
puts stdout "Processing [llength $full_groups_files] group files...\n"

# For each group fasta file, extract just the GPA sequence. We will run the COG detection on the GPA proteins only
set counter			1
set plasmid_counter	0

# set full_groups_files [lrange $full_groups_files 0 100]
foreach group_file $full_groups_files {

	set file_name	[string trim [file tail $group_file]]
	set allProts	[split_genes [string trim [openfile $group_file]]]

	set agProts		{}

	foreach protSeq $allProts {
		set protID	[string trim [string range $protSeq 1 [string first \n $protSeq]]]
		set dbGet	[all_prot_db eval {select is_ag, plasmid from t1 where protID = $protID}]
		set isGeo	[lindex $dbGet 0]
		if {$isGeo == 1} {
			set isPlasmid	[lindex $dbGet 1]
			if {$isPlasmid eq "T"} {
				incr plasmid_counter
			}
			lappend agProts $protSeq
		}
		
	}
	set out		[open $output_dir/$file_name w]
	puts $out	[join $agProts \n]
	close $out

	puts -nonewline "\rDone: $counter / [llength $full_groups_files]"
	flush stdout
	incr counter
}

all_prot_db close