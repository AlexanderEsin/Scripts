#!/usr/local/bin/tclsh
###################################
source  /Users/aesin/Documents/Scripts/General_utils.tcl
package require sqlite3

set direct 			/Users/aesin/Desktop/Geo_again
set final_grp_dir	$direct/Family_mapping/Reduced_groups_fasta
set align_dir		$direct/Group_alignment
file mkdir			$align_dir

set geo_taxid_file	$direct/Genomes/Genome_lists/AG_taxids.txt
set geo_taxids		[split [string trim [openfile $geo_taxid_file]] \n]


## // ##

set final_group_fastas	[glob $final_grp_dir/*fasta]
set test_fasta			[lindex $final_group_fastas 0]


set prot_fastas		[split_genes [string trim [openfile $test_fasta]]]
set ag_prots		{}
foreach prot $prot_seqs {
	set protID		[string range $prot 1 [string first \n $prot]-1]
	set prot_seq	[string range $prot [string first \n $prot]+1 end]

	## Remove all newlines in fasta to get sequence length
	set prot_len	[string length [string map [list \n {}] $prot_seq]]

	## Find whether this protID is from AnoGeo
	set prot_taxid	[string range $protID [string last \. $protID]+1 [string last \_ $protID]-1]
	if {[lsearch $geo_taxids $prot_taxid] != -1} {
		lappend ag_prots $prot
	}
}

















exec clustalo --auto -i $test_fasta -o $align_dir/test_out.fas --outfmt=fasta -v --force -l $align_dir/MSA_log.txt --threads=4