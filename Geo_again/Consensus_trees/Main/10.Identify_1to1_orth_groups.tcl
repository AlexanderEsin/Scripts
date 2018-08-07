#!/usr/local/bin/tclsh

## In this intermediate script we identify family groups
## that represent 1-to-1 orthologous gene families. These
## can then be used to build the supermatrix tree.

source ~/Documents/Scripts/General_utils.tcl
package require sqlite3

set ival			2.0
# set evalue		1e-10
set group_name		[lindex $argv 0]
set evalue			[lindex $argv 1]
set trunc_eval		[string range $evalue 2 end]

## Define the directories and create new output folder
set direct			/users/aesin/Desktop/Geo_again
set consens_dir		$direct/Consensus_groups/$group_name
set prot_db_file	$direct/All_prot_db
set taxid_compare	$direct/Genomes/Genome_lists/Acc_ass_taxid_table.tsv

set group_output	$consens_dir/Family_groups/$trunc_eval
set final_grp_file	$group_output/Final_groups.tsv

## Initialise the log
set log_file		$group_output/Log.txt
set log_out			[open $log_file a]

## Open the protein database
sqlite3 db1 $prot_db_file

## // ##

## First we get a list of the unique taxids for each species
## in our group of taxa. We take a list of the proteomes used
## to assemble our dataset, and convert the acc_ass to taxids
set acc_ass_file	$consens_dir/$group_name\_acc_ass_list.txt
set acc_ass_list	[split [string trim [openfile $acc_ass_file]] \n]
set acc_ass_transl	[split [string trim [openfile $taxid_compare]] \n]

set all_taxids		{}
foreach acc_ass $acc_ass_list {
	set row		[lsearch -inline $acc_ass_transl "$acc_ass\t*"]
	set taxid	[lindex [split $row \t] 1]
	lappend all_taxids $taxid
}

## // ##

## Next we take all the groups derived from MCL and
## reseeded with duplicates/paralogs (Final_groups.tsv)
## and check each remaining group for 1-to-1 orthology

set final_groups	[split [string trim [openfile $final_grp_file]] \n]
set all_taxids_1to1	0
set 1to1_group_list	{}

foreach group $final_groups {
	set group_split		[split $group \t]
	# The first column is the group number
	set protID_list		[lrange $group_split 1 end]

	## 1-to-1 ortholog groups will always have:
	## number of proteins = number of taxids
	if {[llength $protID_list] != [llength $all_taxids]} {
		continue
	}

	## However, num proteins = num taxids does not
	## necessarily imply 1-to-1 orthology; 1 taxid could
	## be mising and there could be 2x proteins from another
	set group_taxids	{}
	foreach protID $protID_list {
		set prot_taxid		[string range $protID [string last \. $protID]+1 [string last \_ $protID]-1]
		lappend group_taxids $prot_taxid
	}
	## Find number of unique taxids per family
	set unique_grp_taxids	[lsort -unique $group_taxids]
	
	## If the number of unique taxids = num expected taxids
	## then we have 1-to-1 orthology
	if {[llength $unique_grp_taxids] == [llength $all_taxids]} {
		incr all_taxids_1to1
		lappend 1to1_group_list	$group
	}	
}

## Write out all the gene families that are 1-to-1 orthologous
set out		[open $group_output/Groups_1to1_orths.tsv w]
puts $out	[join $1to1_group_list \n]
close $out

multiputs stdout $log_out	"\n$all_taxids_1to1 1-to-1 orthologous groups identified with all genomes considered"

## // We get 1143 ortologous gene families (1-to-1) with a 55 genome input // ##
## Now check check whether removing any single taxid significantly increases
## the number of such orthologous groups detected - this could be a sign of 
## poor genome annotation or quality.

proc RemoveOneTaxid {taxid_list group_list previous_max log_chan} {

	set test_table	{"Taxid removed\tNew orth groups"}

	foreach taxid $taxid_list {

		set leftover_taxids	[lremove $taxid_list $taxid]
		set definite_1to1 0
		## One at a time, remove a single gene

		puts "Testing $taxid ..."
		foreach group $group_list {

			## Get a list of the original protIDs in each group
			set group_split		[split $group \t]
			set protID_list		[lrange $group_split 1 end]

			## Using the testing taxid, remove all protIDs in this
			## group that correspond to this taxid
			set prots_to_rem	[lsearch -all -inline $protID_list "*\.$taxid\_*"]

			## Remove each protein carrying that taxid
			set pruned_protIDs	$protID_list	
			foreach prot $prots_to_rem {
				set pruned_protIDs	[lremove $pruned_protIDs $prot]
			}	

			if {[llength $pruned_protIDs] != [llength $leftover_taxids]} {
				continue
			}

			set group_taxids	{}
			foreach protID $pruned_protIDs {
				set prot_taxid		[string range $protID [string last \. $protID]+1 [string last \_ $protID]-1]
				lappend group_taxids $prot_taxid
			}
			set unique_grp_taxids	[lsort -unique $group_taxids]
			
			if {[llength $unique_grp_taxids] == [llength $leftover_taxids]} {
				incr definite_1to1
			}
		}

		set improvement [expr $definite_1to1 - $previous_max]

		lappend test_table	$taxid\t$improvement
	}

	puts $log_chan	[join $test_table \n]
	puts stdout		"\nResults written to log file"
}
multiputs stdout $log_out	"\nWhen removing one taxid at a time we get more orthologous groups:"
RemoveOneTaxid $all_taxids $final_groups $all_taxids_1to1 $log_out

## We find two outlier taxids: 
## 1817674	: would increase the number of orth groups by	: 106
## 495036	: would increase the number of orth groups by	: 49
## All other taxids give much lower values
# puts $log_out "\n"
# puts $log_out [db1 eval {select acc_ass, binomial from t1 where taxid = 1817674 limit 1}]
# # GCF_001624605.1_ASM162460v1 {Geobacillus sp. 8} 3903800

# puts $log_out [db1 eval {select acc_ass, binomial from t1 where taxid = 495036 limit 1}]
# # GCF_000173035.1_ASM17303v1 {Geobacillus sp. G11MC16} 3545187

close $log_out
db1 close

## In both cases the genome size is on the larger end, so this
## discrepancy unlikely to be down to a reduced genome.

## For now we will continue with the full set of genomes,
## if the tree at the end looks questionable we can revisit.




















