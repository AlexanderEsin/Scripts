#!/usr/local/bin/tclsh

## In this script we adjust the ortholog groups produced by MCL
## by reseeding the duplicate proteins, processing true paralogs
## and removing any orthologous proteins / groups that are based
## on plasmid proteins.

source ~/Documents/Scripts/General_utils.tcl
package require sqlite3

set ival			2.0
set evalue		[lindex $argv 0]
# set evalue			"1e-150"
set trunc_eval		[string range $evalue 2 end]

## Define the directories and create new output folder
set direct			/users/aesin/Desktop/Geo_again
set prot_db_file	$direct/All_prot_db
set mcl_group_dir	$direct/RBBH/MCL_groups
set paralogs_dir	$direct/AGeo_inparalogs
set geo_taxid_file	$direct/Genomes/Genome_lists/AG_taxids.txt

set group_output	$direct/Family_groups/$trunc_eval
file mkdir			$group_output

## Open the protein database
sqlite3 db1 $prot_db_file

## Read in the MCL groups file
set mcl_groups_file	$mcl_group_dir/MCL_groups_adj_$ival\_$trunc_eval\.txt
if {[file exists $mcl_groups_file] != 1} {
	error "MCL groups file does not exist!"
}
set mcl_group_list	[split [string trim [openfile $mcl_groups_file]] \n]
set geo_taxids		[split [string trim [openfile $geo_taxid_file]] \n]


## // ##

## Here we first deal with adding duplicates back into
## the orthologous groups identified by RBBH + MCL
set after_dupl_groups	{}
set duplicates_add		{"Group\tAcc_ass\tDuplicated_IDs"}

set group_AG_content	{}
set total_orths			0
set total_AG_content	0
set total_dups_readd	0
set missing_loc			{}
set bad_group			{}

set group_index		1
foreach group $mcl_group_list {

	## Prepare variable to hold new group (with duplicates)
	## and count how many orthologs were in the original group
	set group_dup_incl	{}
	set group_AG_prot	0
	set protID_orths	[split $group \t]
	set total_orths		[expr $total_orths + [llength $protID_orths]]

	puts "Working on adding duplicates ... $group_index // [llength $mcl_group_list]"

	## Check each protID for duplicates and add to the group if necessary
	foreach protID $protID_orths {
		
		## Remove the position from the protID to search
		## for duplicates. Embed * into the variable to help
		## sqlite search.
		if {[string length [string range $protID [string first \_ $protID] [string last \_ $protID]]] < 4} {
			lappend missing_loc $protID
			lappend bad_group	$group_index
			continue
		}

		set trim_position	"[string range $protID 0 [expr [string last \_ $protID] - 1]]_*"
		set all_duplicates	[db1 eval {select protID from t1 where protID glob $trim_position}]

		## If we find multiple entries - i.e. duplicated proteins
		## add the protIDs into the family - a new list
		
		## 1.1.18 - for some reasons some groups still have duplicates - something to do
		## with dup identification during proteome processing (e.g. 702 @ -10 contained both duplicates)
		## so, we check that we don't readd the same protein twice
		foreach duplicate $all_duplicates {
			if {[lsearch $group_dup_incl $duplicate] == -1} {
				lappend group_dup_incl	$duplicate
				set taxid	[string range $duplicate [string last \. $duplicate]+1 [string last \_ $duplicate]-1]
				if {[lsearch $geo_taxids $taxid] != -1} {
					incr group_AG_prot
				}
			}
		}

		## If the length of protIDs is > 1, we have duplicates
		## We create a new table to keep track of those readded
		## and count how many have been added.
		if {[llength $all_duplicates] > 1} {
			set acc_ass				[db1 eval {select acc_ass from t1 where protID = $protID}]
			lappend duplicates_add	$group_index\t$acc_ass\t[join $all_duplicates " "]
			set total_dups_readd	[expr $total_dups_readd + [expr [llength $all_duplicates] - 1]]
		}
	}

	## Add the new group (with duplicates) to the global group list
	## We also add the row/group index - this will now be
	## the gene family / group number
	lappend after_dupl_groups	[join [concat $group_index $group_dup_incl] \t]
	
	## Write the number of AG proteins in each group to a global list
	lappend group_AG_content	[join [concat $group_index $group_AG_prot] \t]
	set total_AG_content		[expr $total_AG_content + $group_AG_prot]

	incr group_index
}

set pruned_groups	$after_dupl_groups
set unique_bad_grps	[lsort -unique $bad_group]
foreach bad_grp $unique_bad_grps {
	set grp				[lsearch -inline $after_dupl_groups "$bad_grp\t*"]
	set pruned_groups	[lremove $pruned_groups $grp]
}
puts "Number of original groups: [llength $after_dupl_groups]"
puts "Number of groups with good protIDs: [llength $pruned_groups]"


## Write out a table listing which duplicates
## were readded to what ortholog groups
set out		[open $group_output/Readded_duplicates.tsv w]
puts $out	[join $duplicates_add \n]
close $out

## Write out the AG content per group at each evalue
set out		[open $group_output/AG_content.tsv w]
puts $out	[join $group_AG_content \n]
close $out

## // ##


## Check what happened to the in-paralogs in the MCL
set in_paralog_data		[string trim [openfile $paralogs_dir/True_paralogs\_$trunc_eval.tsv]]
set in_paralog_list		[lrange [split $in_paralog_data \n] 1 end]

## Make the directory structure to record misclustered in-paralogs
set miscluster_dir		$group_output/Misclustered_inparalogs
set per_para_mis_dir	$miscluster_dir/Per_paralog_pair
file mkdir $miscluster_dir $per_para_mis_dir

## A list of groups (group numbers) that we will mask from
## the final list of gene families. These genes have misclustered
## in-paralogs
set all_masked_groups	{}

set same_group_count	0
set diff_group_count	0
set count	 			1

foreach paralog_pair $in_paralog_list {
	set paralogs		[split $paralog_pair \t]
	set para_A			[lindex $paralogs 0]
	set para_B			[lindex $paralogs 1]

	## Search for paralog A
	set para_A_group_index	[lsearch -all $pruned_groups *\t$para_A*]
	set para_B_group_index	[lsearch -all $pruned_groups *\t$para_B*]

	puts "Determing in-paralog fates ... $count // [llength $in_paralog_list]"
	incr count

	## // ##
	if {[string length $para_A_group_index] == 0 && [string length $para_B_group_index] == 0} {
		error "Should both be missing?"
	}
	if {[string length $para_A_group_index] == 0 || [string length $para_B_group_index] == 0} {
		error "Should one be missing?"
	}

	## We expect in-paralogs to belong to the same group
	## If they do not, we will mask these groups
	if {$para_A_group_index == $para_B_group_index} {
		incr same_group_count
		continue	
	}

	## Both IDs are present but belong to different groups - this is 
	## likely a mis-clustering
	incr diff_group_count
	set pair_missed_dir		$per_para_mis_dir/$diff_group_count
	file mkdir				$pair_missed_dir

	## Set lists to write out group info for the misclustered in-paralogs
	set paralogs_separate	{}
	set para_sep_products	{}

	## Write out the paralog pair that misclustered
	set out			[open $pair_missed_dir/Paralog_pair.tsv w]
	puts $out		[join [concat $para_A $para_B] \n]
	close $out

	## 
	foreach group_index [concat $para_A_group_index $para_B_group_index] {

		## Get the group, add it to a global list for masking
		set group	[lindex $pruned_groups $group_index]
		lappend all_masked_groups $group

		## Get the group number (first column)
		set group_numb	[lindex [split $group \t] 0]
		## Get the protIDs (all other columns)
		set protein_IDs	[lrange [split $group \t] 1 end]

		## For each protID in a group, get the product and make a list
		set group_prods	$group_numb
		foreach protID $protein_IDs	{
			set product	[db1 eval {select product from t1 where protID = $protID};]
			lappend 	group_prods $product
		}

		lappend paralogs_separate	$group
		lappend para_sep_products	[join $group_prods \t]		
	}

	## Write out the groups to which the misclustered in paralogs belong
	set out			[open $pair_missed_dir/Paralog_groups.tsv w]
	puts $out		[join $paralogs_separate \n]
	close $out

	## Write out the products of groups to which the 
	## misclustered in paralogs belong
	set out			[open $pair_missed_dir/Paralog_products.tsv w]
	puts $out		[join $para_sep_products \n]
	close $out
}

## Write out all the masked groups to a combined file
set out			[open $miscluster_dir/all_masked_groups.tsv w]
puts $out		[join $all_masked_groups \n]
close $out

set after_masking_groups		$pruned_groups
foreach masked_group $all_masked_groups {
	set after_masking_groups	[lremove $after_masking_groups $masked_group]	
}

## // ##

## Finally, we remove any groups with fewer than
## 4 members - the minimum required to make a gene tree

set after_pruning_groups	{}
set stub_groups				{}
foreach group $after_masking_groups {
	set protIDs	[lrange [split $group \t] 1 end]
	if {[llength $protIDs] < 4} {
		lappend stub_groups $group
	} else {
		lappend after_pruning_groups $group
	}
}

## Write out a list of final (pruned) groups
set out 		[open $group_output/Final_groups.tsv w]
puts $out		[join $after_pruning_groups \n]
close $out

## Write out a list of the stub groups
set out 		[open $group_output/Stub_groups.tsv w]
puts $out		[join $stub_groups \n]
close $out


# exec clear >@ stdout
set log		[open $group_output/Log.txt w]

multiputs stdout $log "Genes before duplicate addition:\t$total_orths"
multiputs stdout $log "Duplicate genes readded:\t$total_dups_readd"
multiputs stdout $log "Total number of AnoGeo prots:\t$total_AG_content"

multiputs stdout $log "\nBefore removing misclustered groups: [llength $after_dupl_groups] groups"
multiputs stdout $log "After removing groups with bad protIDs: [llength $pruned_groups] groups"
multiputs stdout $log "After removing misclustered groups: [llength $after_masking_groups] groups\n"

multiputs stdout $log "\nGroups with 4 or more protein IDs: [llength $after_pruning_groups] groups"
multiputs stdout $log "Groups with fewer than 4 proteins IDs: [llength $stub_groups] groups\n"


close $log
db1 close







