#!/usr/local/bin/tclsh

## In this script we populate the gene families with 
## the protein sequence. We won't use the full headers
## just the full protID. These are our reference fastas.

## We also create a secondary - almost identical output;
## however we use the taxid only as the header for each
## protein. This works in this case because we know each 
## taxid appears exactly once in each gene family.


source ~/Documents/Scripts/General_utils.tcl
package require sqlite3

set evalue			"1e-50"
set trunc_eval		[string range $evalue 2 end]

## Define the directories and create new output folder
set direct			/Users/aesin/Desktop/Staph
set consens_dir		$direct/Prefilter/Staphylococcus
set prot_db_file	$direct/All_prot_db_new

set group_output	$consens_dir/Family_groups/$trunc_eval
set group_file		$group_output/Groups_1to1_orths.tsv

set ref_out_dir		$group_output/Group_fastas_ref
set tax_out_dir		$group_output/Group_fastas_key
file mkdir			$ref_out_dir $tax_out_dir

## Open the protein database
sqlite3 db1 $prot_db_file

## Read in the 1-to-1 ortholog file
set group_data		[string trim [openfile $group_file]]
set groups			[split $group_data \n]
set counter			1

foreach group $groups {

	puts "Sourcing fasta sequence for group ... $counter / [llength $groups]"

	set group_number	[lindex $group 0]
	set protIDs			[lrange $group 1 end]

	set ref_fasta_l	{}
	set tax_fasta_l	{}

	foreach protID $protIDs {
		## Get the protein sequence
		set prot_AA_seq		[db1 eval {select sequence from t1 where protID = $protID}]
		## Get the taxid
		set taxid			[string range $protID [string last \. $protID]+1 [string last \_ $protID]-1]

		set ref_header		">$protID"
		set tax_header		">$taxid"

		set ref_fasta_seq	[join [concat $ref_header $prot_AA_seq] \n]
		set tax_fasta_seq	[join [concat $tax_header $prot_AA_seq] \n]

		lappend ref_fasta_l	$ref_fasta_seq
		lappend tax_fasta_l	$tax_fasta_seq
	}

	## Write out the fasta files
	set out		[open $ref_out_dir/$group_number\.fasta w]
	puts $out	[join $ref_fasta_l \n]
	close $out

	set out		[open $tax_out_dir/$group_number\.fasta w]
	puts $out	[join $tax_fasta_l \n]
	close $out

	incr counter
}

db1 close