#!/usr/local/bin/tclsh
source ~/Documents/Scripts/General_utils.tcl
package require sqlite3

set master_dir		/Users/aesin/Desktop/Bacillus
set output_dir		$master_dir/Functional_annotation/Core_only_groups
set famGrp_dir		$master_dir/Family_mapping/Reduced_groups_fasta
set protDB_path		$master_dir/ForBac_prot_db


set core_taxid_file		$master_dir/Core_genomes/Genome_lists/coreToKeep.tsv
set core_taxid_tbl		[lrange [split [string trim [openfile $core_taxid_file]] \n] 1 end]
set core_taxids		{}
foreach row $core_taxid_tbl {
	set taxid	[lindex [split $row \t] 2]
	lappend core_taxids $taxid
}

file mkdir $output_dir

# Open database
sqlite3 all_prot_db $protDB_path

# List of reduced fasta groups (input to tree building)
set full_groups_files	[glob $famGrp_dir/*fasta]
puts stdout "Processing [llength $full_groups_files] group files...\n"

# For each group fasta file, extract just the core sequence. We will run the COG detection on the core proteins only
set counter			1

# set full_groups_files [lrange $full_groups_files 0 100]
foreach group_file $full_groups_files {

	set file_name	[string trim [file tail $group_file]]
	set allProts	[split_genes [string trim [openfile $group_file]]]

	set coreProts		{}

	foreach protSeq $allProts {
		set protID	[string trim [string range $protSeq 1 [string first \n $protSeq]]]
		set taxid	[all_prot_db eval {select taxid from t1 where protID = $protID}]
		if {[lsearch $core_taxids $taxid] != -1} {
			lappend coreProts $protSeq
		}
		
	}
	set out		[open $output_dir/$file_name w]
	puts $out	[join $coreProts \n]
	close $out

	puts -nonewline "\rDone: $counter / [llength $full_groups_files]"
	flush stdout
	incr counter
}

all_prot_db close