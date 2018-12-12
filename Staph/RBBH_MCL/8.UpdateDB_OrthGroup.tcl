#!/usr/local/bin/tclsh

source  /Users/aesin/Documents/Scripts/General_utils.tcl
package require sqlite3

set direct 			/Users/aesin/Desktop/Staph

## Define and open all prot DB
set all_db_file		$direct/ForStaph_prot_db
sqlite3 all_prot_db $all_db_file

## Define the directory with final group fasta files
set fmap_dir		$direct/Family_mapping
set grpFasta_dir	$fmap_dir/Reduced_groups_fasta
set grpFasta_list	[lsort -dictionary [glob $grpFasta_dir/*fasta]]

## Add a column to database to hold orthologous group number
puts	"Adding new column to database ..."
all_prot_db eval {ALTER TABLE t1 ADD COLUMN OrthGroup int}
puts -nonewline	"\rAdding new column to database ... done \n"


## For each ortholog group and each protID within it, add the OrthGroup ID to
## the Sqlite table
set counter			1
foreach grpFasta_file $grpFasta_list {
	# Orthologous group number
	set file_name		[file tail $grpFasta_file]
	set grpNumber		[string range $file_name 0 end-6]

	# List of all proteins in group
	set grpProt_list	[split_genes [openfile $grpFasta_file]]
	set grpProt_len		[llength $grpProt_list]

	# Track progress
	puts -nonewline		"\rAdding ortholog group number to database: $counter // [llength $grpFasta_list] ..."
	flush stdout

	# For each protein, get protID and update the Sqlite entry for that protID 
	# with the OrthGroup number
	foreach grpProt $grpProt_list {
		set protID			[string trim [string range $grpProt 1 [string first \n $grpProt]]]
		all_prot_db eval	{UPDATE t1 SET OrthGroup = $grpNumber WHERE protID = $protID}
	}

	incr counter
}

## Make an index on the OrthGroup column
puts "\nMaking index on the OrthGroup column ..."
all_prot_db eval {create index OrthGroup_index on t1 (OrthGroup)}
puts -nonewline "\rMaking index on the OrthGroup column ... done\n"

## Close database
all_prot_db close