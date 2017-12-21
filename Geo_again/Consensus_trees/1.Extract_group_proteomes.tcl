#!/usr/local/bin/tclsh
###################################
source  /Users/aesin/Documents/Scripts/General_utils.tcl
package require sqlite3

set direct 			/Users/aesin/Desktop/Geo_again
set taxdmp_dir		$direct/Genomes/taxdmp
set proteome_dir	$direct/Proteomes
set names_file		$taxdmp_dir/names.dmp
set nodes_file		$taxdmp_dir/nodes.dmp
set accass_tax_file	$direct/Genomes/Genome_lists/Acc_ass_taxid_table.tsv

## // ##
## Desired classification for which to filder genomes
set name			"Bacillales"
puts stdout			"Collecting all proteomes under the taxonomic group: $name ..."
set out_dir			$direct/Consensus_groups
file mkdir			$out_dir/$name

## // ##
## Identify the node corresponding to desired name
set names_data		[split [string trim [openfile $names_file]] \n]
set corr_lines		[lsearch -all -inline $names_data "*\t|\t$name\t|\t*"]

## Make sure we only find one entry
if {[llength $corr_lines] == 1} {
	set name_entry	[lindex $corr_lines 0]
} else {
	error "We find multiple entries for the name: $name"
}

## Process string and extract the node number
set entry_list		[split $name_entry \|]
set node_value		[string trim [lindex $entry_list 0]]
puts stdout			"Identified the taxonomic node value for $name: $node_value"


## // ##
set nodes_data		[split [string trim [openfile $nodes_file]] \n]
set accass_tax_tbl	[split [string trim [openfile $accass_tax_file]] \n]

puts stdout			"Identifying whether genome is part of $name ..."
set group_list		{}
set counter			1
foreach entry $accass_tax_tbl {

	set acc_ass		[lindex [split $entry \t] 0]
	set base_taxid	[lindex [split $entry \t] 1]
	set child_taxid	$base_taxid

	while 1 {
		set entry		[lsearch -inline $nodes_data "$child_taxid\t*"]
		set parent_node	[string trim [lindex [split $entry \|] 1]]

		if {$parent_node == $node_value} {
			lappend group_list	$acc_ass
			break
		} elseif {$parent_node == 1} {
			break
		} else {
			set child_taxid $parent_node
			continue
		}
	}

	puts -nonewline "\tTaxonomic tree scan: $counter // [llength $accass_tax_tbl]\r"
	flush stdout
	incr counter
}

puts stdout			"\nThere are [llength $group_list] genomes in $name ..."
set out				[open $out_dir/$name/$name\_acc_ass_list.txt w]
puts $out 			[join $group_list \n]
close $out

## // ##
## Transfer clean and duplicate fasta files
set group_clean_fasta_dir	$out_dir/$name/Proteome_clean_fastas
set group_dup_fasta_dir		$out_dir/$name/Proteome_dup_fastas
file mkdir $group_clean_fasta_dir $group_dup_fasta_dir

puts stdout 		"Copying proteome fasta files for $name genomes ..."
foreach acc_ass $group_list {
	set fasta_file_name		"$acc_ass\_proteome.fasta"
	file copy -force $proteome_dir/Proteome_clean_fastas/$fasta_file_name $group_clean_fasta_dir
	file copy -force $proteome_dir/Proteome_dup_fastas/$fasta_file_name $group_dup_fasta_dir
}
puts stdout			"All done."






