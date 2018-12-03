#!/usr/local/bin/tclsh
###################################
source  /Users/aesin/Documents/Scripts/General_utils.tcl
package require sqlite3
package require struct::set

set ival			2.0
set evalue_list		[list "-10" "-50" "-100" "-150"]

set direct 			/Users/aesin/Desktop/Bacillus

set all_db_file		$direct/ForBac_prot_db
set core_taxid_file	$direct/Bac_genomes/Genome_lists/coreToKeep.tsv

set group_dir		$direct/Family_groups
set fmap_dir		$direct/Family_mapping

set fmap_db_dir		$fmap_dir/MCL_groups_DBs
set fmap_fmap_dir	$fmap_dir/Fmaps

set out_fasta_dir	$fmap_dir/Reduced_groups_fasta
set out_list		$fmap_dir/Reduced_groups.tsv

## Make necessary directories
file mkdir $fmap_db_dir $fmap_fmap_dir $out_fasta_dir

sqlite3 all_prot_db $all_db_file

set min_gene_context 200


###################################
### Make DBs of MCL groupings at different evalues

foreach evalue $evalue_list {
	set mcl_group_file	$group_dir/$evalue/Final_groups.tsv
	set mcl_group_db	$fmap_db_dir/MCL_group_$ival\_$evalue\_DB

	if {[file exists $mcl_group_db] == 1} {
		puts stdout "DB at evalue : $evalue already exists ... skipping"
		continue
	}

	set mcl_groups		[split [string trim [openfile $mcl_group_file]] \n]
	set num_mcl_groups	[llength $mcl_groups]

	## Prep the DB
	sqlite3 db1 $mcl_group_db
	db1 eval {PRAGMA main.page_size = 4096}
	db1 eval {PRAGMA main.cache_size = 10000}
	db1 eval {PRAGMA main.locking_mode = EXCLUSIVE}
	db1 eval {PRAGMA main.synchronous = NORMAL}
	db1 eval {PRAGMA main.journal_mode = WAL}
	## Two columns, group number and orthologs
	db1 eval {create table t1(group_num int, orthologs text)}

	set group_counter 1
	foreach group $mcl_groups {
		puts -nonewline stdout "Adding orthologs to DB -- evalue : $evalue -- group: $group_counter // $num_mcl_groups\r"
		flush stdout
		set group_split	[split [string trim $group] \t]
		
		## The first element is the group number, the rest orthologs
		set group_num	[lindex $group_split 0]
		set orthologs	[lrange $group_split 1 end]

		foreach orth $orthologs {
			db1 eval {insert into t1 values($group_num,$orth)}
		}
		incr group_counter
	}

	puts stdout "Making index on DB -- evalue : $evalue"
	db1 eval {CREATE INDEX gene_index ON t1 (orthologs)}
	db1 eval {CREATE INDEX grp_index ON t1 (group_num)}

	puts stdout "DB at evalue : $evalue ... DONE\n"
	db1 close
}

###################################
### Link clusters across different evalues

set master_counter	0
set rev_eval_list	[lreverse $evalue_list]
set starting_eval	[lindex $rev_eval_list 0]

set core_taxid_tbl		[split [string trim [openfile $core_taxid_file]] \n]
set core_taxids		{}
foreach row $core_taxid_tbl {
	set taxid	[lindex [split $row \t] 2]
	lappend core_taxids $taxid
}

while {[llength $rev_eval_list] >= 2} {
	if {$master_counter != 0} {
		set starting_eval	[lindex $rev_eval_list 1]
	}

	## At every iteration, remove the highest stringency cut-off and construct
	while {[lindex $rev_eval_list 0] != $starting_eval} {
		set rev_eval_list [lrange $rev_eval_list 1 end]
	}
	set eval_number	[llength $rev_eval_list]
	

	set dir_counter	0
	set conserved_families	{}
	set error_nos	{}
	set family_map	{}


	if {$eval_number >=2} {
		puts "\n\nE-value cutoffs being tested: $rev_eval_list"
	} else {
		puts "\n\nLowest_stringency reference file being written..."
	}

	### When the e-value cutoff can be compared to lower stringencies create a fmaily map ###
	while {$dir_counter <= [expr $eval_number - 2]} {

		set dir_0	[lindex $rev_eval_list $dir_counter]
		set dir_1	[lindex $rev_eval_list [expr $dir_counter + 1]]

		# set counter 0

		## Open the high e-val (e.g. -100) clustering as a database
		sqlite3 lower_db $fmap_db_dir/MCL_group_2.0_$dir_1\_DB

		## Open the lower e-val (e.g. -150) clustering as the query, and split by gene family
		set groups_file			$group_dir/$dir_0/Final_groups.tsv
		set gene_families_0		[split [string trim [openfile $groups_file]] \n]

		# # If this is the first comparison ... #
		if {$dir_counter == 0} {
			set this_round	$gene_families_0

		} else {
			set this_round			{}		
			foreach group_num $conserved_families {
				set family [lsearch -inline -glob $gene_families_0 "$group_num\t*"]
				lappend this_round $family
			}
		}

		## Set the temporary conserved families for this iteration of the comparison
		set conserved_families	{}
		set conserved_temp		{}
		set family_map_temp		{}


		puts "\tNow scanning @ $dir_1"
		foreach family $this_round {

			# At the higher stringency level (e.g. -150) for each group, get all the protIDs #
			set group_num	[lindex [split [string trim $family] \t] 0]
			set prot_IDs	[lrange [split [string trim $family] \t] 1 end]

			# Within those find all the core protIDs #
			set core_IDs		{}
			foreach prot_ID $prot_IDs {
				set taxid	[string range $prot_ID [string last \. $prot_ID]+1 [string last \_ $prot_ID]-1]
				if {[lsearch $core_taxids $taxid] != -1} {
					lappend core_IDs $prot_ID
				}
			}
			
			if {[llength $core_IDs] == 0} {
				error "We expect to find some core_IDs at lower stringency. Cannot find core prot IDs @ $dir_0 in group: $group_num"
			}
		
			# Foreach core protID, identify the corresponding group in the lower stringency level (e.g. -100)
			set index_find [list]
			foreach core_ID $core_IDs {
				# Query the lower stringency
				set lower_index	[lower_db eval {select group_num from t1 where orthologs = $core_ID}]

				# A core_ID might not be found at a lower stringency because it was part of
				# misclustered InParalog groups that have been masked out. Alternatively,
				# it could be clustered in a stub group (<4 protIDs) that was filtered out.
				if {[string length $lower_index] == 0} {
					# puts stdout "core_ID $core_ID @ $dir_0 in group: $group_num not found in any group @ $dir_1"
					continue
				}

				# List of all the corresponding group numss from the lower stringency level for each core protID 
				lappend index_find $lower_index
			}

			# Unique group numbers to which these core_IDs belong at a lower stringency level
			set unique_indices	[lsort -unique $index_find]
			
			# If all of the group nos map to a single group (i.e. the core prots have not split into two groups at the lower stringency level), we add the index (lower stringency group no) to the conserved temp list. A map is also constructed (original higher stringency -> new index lower stringency), and added to the family_map_temp variable which is global for two stringency comparisons
			if {[llength $unique_indices] == 1} {
				lappend conserved_temp	$unique_indices
				lappend family_map_temp	[join [list $group_num $unique_indices] \t]
			}
		}

		# Once all the indices for the conserved groups have been collected - we search for duplicates. Duplicates (e.g. two -150 families mapped exclusively to one -100 family) suggest coalescence #
		# Any entries that have coalesced into one group at the lower stringency are removed from the lists #
		set coalescent_families [dups $conserved_temp]
		foreach family $conserved_temp {
			if {[lsearch $coalescent_families $family] == -1} {
				lappend conserved_families $family
			}
		}

		if {$dir_counter == 0} {
			foreach family $family_map_temp {
				set destination_family [lindex [split $family \t] 1]
				if {[lsearch $coalescent_families $destination_family] == -1} {
					lappend family_map $family
				}
			}
		} else {
			set new_family_map		{}
			foreach family $family_map_temp {
				set present_family		[lindex $family 0]
				set destination_family	[lindex $family 1]
				## If the destination family is NOT coalescent then add to new family map
				if {[lsearch $coalescent_families $destination_family] == -1} {
					set existing_entry	[split [lsearch -inline -glob $family_map "*\t$present_family"] \t]
					set new_entry		[join [concat $existing_entry $destination_family] \t]
					lappend new_family_map $new_entry
				}
			}
			set family_map $new_family_map
		}
	
		lower_db close
		incr dir_counter
	}

	## Output the family maps at every e-value cutoff except the lowest one
	if {$eval_number >= 2} {
		set fam_map_out	[join $family_map \n]
		set out			[open $fmap_fmap_dir/Fmap_$starting_eval\.tsv w]
		puts $out 		$fam_map_out
		close $out

		incr master_counter
	}

	### Make the reference list at the lowest e-value stringency ###
	if {$eval_number == 1} {

		set reference_index {}

		set lowest_dir [lindex $evalue_list 0]
		set groups_file			$group_dir/$lowest_dir/Final_groups.tsv
		set gene_families_low	[split [string trim [openfile $groups_file]] \n]

		foreach family $gene_families_low {

			set group_num	[lindex [split [string trim $family] \t] 0]
			set prot_IDs	[lrange [split [string trim $family] \t] 1 end]

			set core_IDs		{}
			foreach prot_ID $prot_IDs {
				set taxid	[string range $prot_ID [string last \. $prot_ID]+1 [string last \_ $prot_ID]-1]
				if {[lsearch $core_taxids $taxid] != -1} {
					lappend core_IDs $prot_ID
				}
			}
			
			if {[llength $core_IDs] == 0} {
				error "We expect to find some core_IDs in each family. Cannot find core prot IDs @ $dir_0 in group: $group_num"
			}

			lappend family_map $group_num
		}

		set fam_map_out	[join $family_map \n]
		set out			[open $fmap_fmap_dir/Fmap_$starting_eval\.tsv w]
		puts $out 		$fam_map_out
		close $out

		incr master_counter
	}		
}



###################################
### Pick best cluster for each basal group

set rev_eval_list	[lreverse $evalue_list]

set base_evalue		[lindex $rev_eval_list end]
set base_core_content	[split [string trim [openfile $group_dir/$base_evalue/Core_prot_content.tsv]] \n]
set base_group_list	[split [string trim [openfile $group_dir/$base_evalue/Final_groups.tsv]] \n]
set base_group_tot	[llength $base_group_list]

set base_fams_done		{}
set done_counter		1
set family_group_table	[list [join [list "Base group" "Base protID num" "Evalue sampled" "At evalue group" "At evalue protID num" "Base core number" "Cluster core number"] \t]]

foreach eval $rev_eval_list {
	
	set fmap_data		[split [string trim [openfile $fmap_fmap_dir/Fmap_$eval\.tsv]] \n]
	set test_fmap_len	[split [lindex $fmap_data 0] \t]
	set fmap_index		[expr [llength $rev_eval_list] - [llength $test_fmap_len]]
	set test_evals		[lrange $rev_eval_list $fmap_index end]
		
	foreach fmap $fmap_data {

		set group_clusters	[split $fmap \t]
		## The families are eventually named based on the group number at the most basal (e.g. -10) level ##
		set base_group_num	[lindex $group_clusters end]

		## If the base group has already been done (e.g. at a higher fmap)
		## skip and continue
		if {[lsearch -exact $base_fams_done $base_group_num] != -1} {
			continue
		} else {

			## A cluster picked at non-basal evalue must contain all the basal core sequences
			## Count how many core sequences are represented at basal evalue
			set min_core_required	[lindex [split [lsearch -inline -glob $base_core_content "$base_group_num\t*"] \t] 1]

			set cluster_counter 0
			## The fmaps run lowest evalue up (e.g. -150 -100 -50 -10). We query each starting from the most stringent
			foreach cluster $group_clusters {

				set test_evalue		[lindex $test_evals $cluster_counter]

				## If this is the basal evalue (e.g. -10) we take the family as is
				## with no further checks
				if {$test_evalue == $base_evalue} {
					## Get the base number of proteins
					set base_group		[lsearch -inline -glob $base_group_list "$base_group_num\t*"]
					set base_protIDs	[lrange [split $base_group \t] 1 end]
					set base_prot_num	[llength $base_protIDs]

					## Find and collate the fasta sequences for the picked group
					set group_fasta		{}
					foreach base_protID $base_protIDs {
						set prot_data	[all_prot_db eval {select protID, sequence from t1 where protID = $base_protID}]
						set fasta_prot	">[join $prot_data \n]"
						lappend group_fasta $fasta_prot
					}

					## Write out the final group fasta
					set out		[open $out_fasta_dir/$base_group_num\.fasta w]
					puts $out	[join $group_fasta \n]
					close $out

					## Append global lists
					lappend base_fams_done		$base_group_num
					lappend family_group_table	[join [list $base_group_num $base_prot_num $test_evalue $cluster $base_prot_num $min_core_required $min_core_required] \t]

					puts stdout	"DONE: $done_counter /// $base_group_tot\t[join [list $base_group_num $base_prot_num $test_evalue $cluster $base_prot_num $min_core_required] \t]"
					incr done_counter
					break
				}


				## How many core in this cluster
				set test_core_content	[split [string trim [openfile $group_dir/$test_evalue/Core_prot_content.tsv]] \n]
				set clust_core_cont	[lindex [split [lsearch -inline -glob $test_core_content "$cluster\t*"] \t] 1]

				## If there are not enough core genes in the mapped cluster at this evalue, continue to the next one
				if {$clust_core_cont < $min_core_required} {
					incr cluster_counter
					continue
				}

				## Otherwise check how many genes there are in this cluster
				sqlite3 test_db $fmap_db_dir/MCL_group_$ival\_$test_evalue\_DB
				set clust_prot_IDs	[test_db eval {select orthologs from t1 where group_num = $cluster}]
				test_db close

				set clust_prot_num	[llength $clust_prot_IDs]

				## Continue if the number of genes is less than the minimum wanted context
				if {$clust_prot_num < $min_gene_context} {
					incr cluster_counter
					continue
				}

				set base_protIDs	[lsearch -inline -glob $base_group_list "$base_group_num\t*"]
				set base_prot_num	[llength [lrange [split $base_protIDs \t] 1 end]]

				## Find and collate the fasta sequences for the picked group
				set group_fasta		{}
				foreach clust_protID $clust_prot_IDs {
					set prot_data	[all_prot_db eval {select protID, sequence from t1 where protID = $clust_protID}]
					set fasta_prot	">[join $prot_data \n]"
					lappend group_fasta $fasta_prot
				}

				## Write out the final group fasta
				set out		[open $out_fasta_dir/$base_group_num\.fasta w]
				puts $out	[join $group_fasta \n]
				close $out

				lappend base_fams_done		$base_group_num
				lappend family_group_table	[join [list $base_group_num $base_prot_num $test_evalue $cluster $clust_prot_num $min_core_required $clust_core_cont] \t]
				puts stdout	"DONE: $done_counter /// $base_group_tot\t[join [list $base_group_num $base_prot_num $test_evalue $cluster $clust_prot_num $clust_core_cont] \t]"
				incr done_counter
				break
			}
		}
	}
}

set out		[open $out_list w]
puts $out	[join $family_group_table \n]
close $out

all_prot_db close

# ###########################################################################

# puts "ALL DONE!"