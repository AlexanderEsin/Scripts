#!/usr/local/bin/tclsh
###########################################################################
## Changelog ##
# 23 November 2015 - changed "Error 2" to "Exit 2" so that the script actually breaks at those points #
# 14 January 2016 - added an option to pick the most stringent evalue at which all the Geobacillus sequences are present #


###########################################################################
source  /Users/aesin/Documents/Scripts/General_utils.tcl
package require sqlite3

set ival			2.0
set evalue_list		[list "-10" "-50" "-100" "-150"]

set direct 			/Users/aesin/Desktop/Geo_again
set group_dir		$direct/Family_groups
set fmap_dir		$direct/Family_mapping

set fmap_db_dir		$fmap_dir/MCL_groups_DBs
set fmap_fmap_dir	$fmap_dir/Fmaps
set output_file		$fmap_dir/Reduced_gene_families

file mkdir $fmap_db_dir $fmap_fmap_dir


# set make_db_switch 0
# set make_fmap_switch 0
# set final_groups_switch 1
# ###
# set minimum_gene_context 50
# set minimum_geo_number "all"; # Acceptable values: integer >= 1 or "all" #
# ###


#########################

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
		puts stdout "Adding orthologs to DB -- evalue : $evalue -- group: $group_counter // $num_mcl_groups"
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

	puts stdout "DB at evalue : $evalue ... DONE\n"
	db1 close
}

# #########################


# set master_counter	0
# set rev_eval_list	[lreverse $evalue_list]
# set starting_eval	[lindex $evalue_list 0]


# while {[llength $rev_eval_list] >= 2} {
# 	if {$master_counter != 0} {
# 		set starting_eval	[lindex $rev_eval_list 1]
# 	}

# 	## At every iteration, remove the highest stringency cut-off and construct
# 	while {[lindex $rev_eval_list 0] != $starting_eval} {
# 		set rev_eval_list [lrange $rev_eval_list 1 end]
# 	}
# 	set eval_number	[llength $rev_eval_list]
	

# 	set dir_counter	0
# 	set conserved_families	{}
# 	set error_nos	{}
# 	set family_map	{}


# 	if {$eval_number >=2} {
# 		puts "\n\nE-value cutoffs being tested: $rev_eval_list"
# 	} else {
# 		puts "\n\nLowest_stringency reference file being written..."
# 	}


# 	### When the e-value cutoff can be compared to lower stringencies create a fmaily map ###
# 	while {$dir_counter <= [expr $eval_number - 2]} {

# 		set dir_0	[lindex $rev_eval_list $dir_counter]
# 		set dir_1	[lindex $rev_eval_list [expr $dir_counter + 1]]

# 		# set counter 0

# 		## Open the high e-val (e.g. -100) clustering as a database
# 		sqlite3 db1 $fmap_db_dir/MCL_group_2.0_$dir_1\_DB

# 		## Open the lower e-val (e.g. -150) clustering as the query, and split by gene family
# 		cd $direct/$ival/$dir_0
# 		openfile Group_list_$ival\.tsv
# 		set gene_families_0 [split [string trim $data] \n]

# 		### Set the temporary conserved families for this iteration of the comparison ###
# 		set conserved_temp {}
# 		set error_temp {}
# 		set family_map_temp {}

# 		# If this is the first comparison ... #
# 		if {$dir_counter == 0} {
# 			puts "\tNow scanning $dir_1"
# 			foreach family $gene_families_0 {
# 				set map_entry {}
# 				# At the higher stringency level (e.g. -150) for each gene family, get all the gene ids #
# 				set genes [lrange [split [string trim $family] \t] 1 end]
# 				# Within those find all the geobacillus genes #
# 				set geo_genes [lsearch -all -glob -inline $genes *\(Geobac\)]
# 				# If there are no geobacillus genes, skip to the next gene family #
# 				if {[llength $geo_genes] == 0} {
# 					continue
# 				}

# 				set group_number [string trim [string range $family 0 [string first \t $family]]]

# 				set index_find {}
# 				# Foreach geobacillus gene, identify the corresponding gene family in the lower stringency level (e.g. -100) #
# 				foreach gene $geo_genes {
# 					set index ""
# 					# Query the lower stringency directory #
# 					db1 eval {SELECT group_no FROM t1 WHERE genes = $gene} {
# 						set index "$group_no"
# 					}
# 					if {$index == ""} {
# 						puts "Group number not found ...."
# 						exit 2
# 					}
# 					# List of all the corresponding group nos from the lower stringency level for each geobacillus gene #
# 					lappend index_find $index
# 				}

# 				# Identify the corresponding group nos without repeats #
# 				set unique_indices [lsort -unique $index_find]
				
# 				# If all of the group nos map to a single gene family (i.e. the geobacilli have not split into two groups at the lower stringency level), we add the index (lower stringency group no) to the conserved temp list. A map is also constructed (original higher stringency -> new index lower stringency), and added to the family_map_temp variable which is global for two stringency comparisons #
# 				if {[llength $unique_indices] == 1} {
# 					lappend conserved_temp $unique_indices

# 					lappend map_entry $group_number $unique_indices
# 					set map_entry [join $map_entry \t]
# 					lappend family_map_temp $map_entry
# 				}
# 			}

# 			# Once all the indices for the conserved gene families have been collected - we search for duplicates. Duplicates (e.g. two -150 families mapped exclusively to one -100 family) suggest coalescence #
# 			# Any entries that have coalesced into one group at the lower stringency are removed from the lists #
# 			set coalescent_families [dups $conserved_temp]
# 			foreach family $conserved_temp {
# 				if {[lsearch $coalescent_families $family] == -1} {
# 					lappend conserved_families $family
# 				}
# 			}
# 			foreach family $family_map_temp {
# 				set destination_family [lindex $family 1]
# 				if {[lsearch $coalescent_families $destination_family] == -1} {
# 					lappend family_map $family
# 				}
# 			}
# 		} else {
# 			puts "\tNow scanning $dir_1"
# 			set next_round {}
# 			set new_family_map {}
# 			foreach id $conserved_families {
# 				set family [lsearch -all -inline -glob $gene_families_0 "$id\t*"]
# 				lappend next_round $family
# 			}

# 			set conserved_families {}
# 			foreach family $next_round {
# 				set family [join $family]
# 				set map_entry {}
# 				set genes [lrange [split [string trim $family] \t] 1 end]
# 				set geo_genes [lsearch -all -glob -inline $genes *\(Geobac\)]

# 				if {[llength $geo_genes] == 0} {
# 					continue
# 				}

# 				set group_number [string trim [string range $family 0 [string first \t $family]]]

# 				set index_find {}
# 				foreach gene $geo_genes {
# 					set index ""
# 					db1 eval {SELECT group_no FROM t1 WHERE genes = $gene} {
# 						set index "$group_no"
# 					}
# 					if {$index == ""} {
# 						puts "Group number not found ...."
# 						exit 2
# 					}

# 					lappend index_find $index
# 				}

# 				set unique_indices [lsort -unique $index_find]
# 				if {[llength $unique_indices] == 1} {
# 					lappend conserved_temp $unique_indices

# 					lappend map_entry $group_number $unique_indices
# 					set map_entry [join $map_entry \t]
# 					lappend family_map_temp $map_entry
# 				}
# 			}
# 			set coalescent_families [dups $conserved_temp]
# 			foreach family $conserved_temp {
# 				if {[lsearch $coalescent_families $family] == -1} {
# 					lappend conserved_families $family
# 				}
# 			}
# 			foreach family $family_map_temp {
# 				set present_family [lindex $family 0]
# 				set destination_family [lindex $family 1]
# 				if {[lsearch $coalescent_families $destination_family] == -1} {
# 					set existing_entry [lsearch -all -inline -glob $family_map "*\t$present_family"]
# 					set existing_entry [join $existing_entry \t]
# 					set new_entry "$existing_entry\t$destination_family"
# 					lappend new_family_map $new_entry
# 				}
# 			}
# 			set family_map $new_family_map
# 		}
# 		#puts $dir_counter
# 		incr dir_counter
# 	}

# 	### Output the family maps at every e-value cutoff except the lowest one ###
# 	if {$eval_number >= 2} {
# 		set fam_map_str [join $family_map \n]
# 		set out [open $fmap_path/Fmaps/Fmap_$starting_eval\.tsv w]
# 		puts $out $fam_map_str
# 		close $out

# 		puts "Master counter: [expr $master_counter + 1]"
# 		incr master_counter
# 	}

# 	### Make the reference list at the lowest e-value stringency ###
# 	if {$eval_number == 1} {

# 		set family_map {}

# 		set dir_0 [lindex $eval_dirs 0]
# 		cd $direct/$ival/$dir_0
# 		openfile Group_list_$ival\.tsv
# 		set gene_families_0 [split [string trim $data] \n]

# 		foreach family $gene_families_0 {

# 			set genes [lrange [split [string trim $family] \t] 1 end]
# 			set geo_genes [lsearch -all -glob -inline $genes *\(Geobac\)]
# 			if {[llength $geo_genes] == 0} {
# 				continue
# 			}

# 			set group_number [string trim [string range $family 0 [string first \t $family]]]
# 			lappend family_map $group_number
# 		}

# 		set fam_map_str [join $family_map \n]
# 		set out [open $fmap_path/Fmaps/Fmap_$starting_eval\.tsv w]
# 		puts $out $fam_map_str
# 		close $out

# 		puts "Master counter: [expr $master_counter + 1]"
# 		incr master_counter
# 	}

			
# }



# # if {$make_fmap_switch == 1} {
# # 	file mkdir $fmap_path/Fmaps

# # 	set master_counter 0
# # 	set starting_eval [lindex $eval_dirs 0]
# # 	while {[llength $eval_dirs] >= 2} {
# # 		if {$master_counter != 0} {
# # 			set starting_eval [lindex $eval_dirs 1]
# # 		}
		
# # 		### At every iteration, remove the highest stringency cut-off and construct ###
# # 		while {[lindex $eval_dirs 0] != $starting_eval} {
# # 			set eval_dirs [lrange $eval_dirs 1 end]
# # 		}
# # 		set eval_number [llength $eval_dirs]
		

# # 		set dir_counter 0
# # 		set conserved_families {}
# # 		set error_nos {}

# # 		set family_map {}

# # 		if {$eval_number >=2} {
# # 			puts "\n\nE-value cutoffs being tested: $eval_dirs"
# # 		} else {
# # 			puts "\n\nLowest_stringency reference file being written..."
# # 		}

# # 		### When the e-value cutoff can be compared to lower stringencies create a fmaily map ###
# # 		while {$dir_counter <= [expr $eval_number - 2]} {

# # 			set dir_0 [lindex $eval_dirs $dir_counter]
# # 			set dir_1 [lindex $eval_dirs [expr $dir_counter + 1]]

# # 			# set counter 0

# # 			### Open the high e-val (e.g. -100) clustering as a database ###
# # 			cd $direct/$ival/$dir_1
# # 			sqlite3 db1 Groups_$dir_1\_db
# # 			### Open the lower e-val (e.g. -150) clustering as the query, and split by gene family ###
# # 			cd $direct/$ival/$dir_0
# # 			openfile Group_list_$ival\.tsv
# # 			set gene_families_0 [split [string trim $data] \n]

# # 			### Set the temporary conserved families for this iteration of the comparison ###
# # 			set conserved_temp {}
# # 			set error_temp {}
# # 			set family_map_temp {}

# # 			# If this is the first comparison ... #
# # 			if {$dir_counter == 0} {
# # 				puts "\tNow scanning $dir_1"
# # 				foreach family $gene_families_0 {
# # 					set map_entry {}
# # 					# At the higher stringency level (e.g. -150) for each gene family, get all the gene ids #
# # 					set genes [lrange [split [string trim $family] \t] 1 end]
# # 					# Within those find all the geobacillus genes #
# # 					set geo_genes [lsearch -all -glob -inline $genes *\(Geobac\)]
# # 					# If there are no geobacillus genes, skip to the next gene family #
# # 					if {[llength $geo_genes] == 0} {
# # 						continue
# # 					}

# # 					set group_number [string trim [string range $family 0 [string first \t $family]]]

# # 					set index_find {}
# # 					# Foreach geobacillus gene, identify the corresponding gene family in the lower stringency level (e.g. -100) #
# # 					foreach gene $geo_genes {
# # 						set index ""
# # 						# Query the lower stringency directory #
# # 						db1 eval {SELECT group_no FROM t1 WHERE genes = $gene} {
# # 							set index "$group_no"
# # 						}
# # 						if {$index == ""} {
# # 							puts "Group number not found ...."
# # 							exit 2
# # 						}
# # 						# List of all the corresponding group nos from the lower stringency level for each geobacillus gene #
# # 						lappend index_find $index
# # 					}

# # 					# Identify the corresponding group nos without repeats #
# # 					set unique_indices [lsort -unique $index_find]
					
# # 					# If all of the group nos map to a single gene family (i.e. the geobacilli have not split into two groups at the lower stringency level), we add the index (lower stringency group no) to the conserved temp list. A map is also constructed (original higher stringency -> new index lower stringency), and added to the family_map_temp variable which is global for two stringency comparisons #
# # 					if {[llength $unique_indices] == 1} {
# # 						lappend conserved_temp $unique_indices

# # 						lappend map_entry $group_number $unique_indices
# # 						set map_entry [join $map_entry \t]
# # 						lappend family_map_temp $map_entry
# # 					}
# # 				}

# # 				# Once all the indices for the conserved gene families have been collected - we search for duplicates. Duplicates (e.g. two -150 families mapped exclusively to one -100 family) suggest coalescence #
# # 				# Any entries that have coalesced into one group at the lower stringency are removed from the lists #
# # 				set coalescent_families [dups $conserved_temp]
# # 				foreach family $conserved_temp {
# # 					if {[lsearch $coalescent_families $family] == -1} {
# # 						lappend conserved_families $family
# # 					}
# # 				}
# # 				foreach family $family_map_temp {
# # 					set destination_family [lindex $family 1]
# # 					if {[lsearch $coalescent_families $destination_family] == -1} {
# # 						lappend family_map $family
# # 					}
# # 				}
# # 			} else {
# # 				puts "\tNow scanning $dir_1"
# # 				set next_round {}
# # 				set new_family_map {}
# # 				foreach id $conserved_families {
# # 					set family [lsearch -all -inline -glob $gene_families_0 "$id\t*"]
# # 					lappend next_round $family
# # 				}

# # 				set conserved_families {}
# # 				foreach family $next_round {
# # 					set family [join $family]
# # 					set map_entry {}
# # 					set genes [lrange [split [string trim $family] \t] 1 end]
# # 					set geo_genes [lsearch -all -glob -inline $genes *\(Geobac\)]

# # 					if {[llength $geo_genes] == 0} {
# # 						continue
# # 					}

# # 					set group_number [string trim [string range $family 0 [string first \t $family]]]

# # 					set index_find {}
# # 					foreach gene $geo_genes {
# # 						set index ""
# # 						db1 eval {SELECT group_no FROM t1 WHERE genes = $gene} {
# # 							set index "$group_no"
# # 						}
# # 						if {$index == ""} {
# # 							puts "Group number not found ...."
# # 							exit 2
# # 						}

# # 						lappend index_find $index
# # 					}

# # 					set unique_indices [lsort -unique $index_find]
# # 					if {[llength $unique_indices] == 1} {
# # 						lappend conserved_temp $unique_indices

# # 						lappend map_entry $group_number $unique_indices
# # 						set map_entry [join $map_entry \t]
# # 						lappend family_map_temp $map_entry
# # 					}
# # 				}
# # 				set coalescent_families [dups $conserved_temp]
# # 				foreach family $conserved_temp {
# # 					if {[lsearch $coalescent_families $family] == -1} {
# # 						lappend conserved_families $family
# # 					}
# # 				}
# # 				foreach family $family_map_temp {
# # 					set present_family [lindex $family 0]
# # 					set destination_family [lindex $family 1]
# # 					if {[lsearch $coalescent_families $destination_family] == -1} {
# # 						set existing_entry [lsearch -all -inline -glob $family_map "*\t$present_family"]
# # 						set existing_entry [join $existing_entry \t]
# # 						set new_entry "$existing_entry\t$destination_family"
# # 						lappend new_family_map $new_entry
# # 					}
# # 				}
# # 				set family_map $new_family_map
# # 			}
# # 			#puts $dir_counter
# # 			incr dir_counter
# # 		}

# # 		### Output the family maps at every e-value cutoff except the lowest one ###
# # 		if {$eval_number >= 2} {
# # 			set fam_map_str [join $family_map \n]
# # 			set out [open $fmap_path/Fmaps/Fmap_$starting_eval\.tsv w]
# # 			puts $out $fam_map_str
# # 			close $out

# # 			puts "Master counter: [expr $master_counter + 1]"
# # 			incr master_counter
# # 		}

# # 		### Make the reference list at the lowest e-value stringency ###
# # 		if {$eval_number == 1} {

# # 			set family_map {}

# # 			set dir_0 [lindex $eval_dirs 0]
# # 			cd $direct/$ival/$dir_0
# # 			openfile Group_list_$ival\.tsv
# # 			set gene_families_0 [split [string trim $data] \n]

# # 			foreach family $gene_families_0 {

# # 				set genes [lrange [split [string trim $family] \t] 1 end]
# # 				set geo_genes [lsearch -all -glob -inline $genes *\(Geobac\)]
# # 				if {[llength $geo_genes] == 0} {
# # 					continue
# # 				}

# # 				set group_number [string trim [string range $family 0 [string first \t $family]]]
# # 				lappend family_map $group_number
# # 			}

# # 			set fam_map_str [join $family_map \n]
# # 			set out [open $fmap_path/Fmaps/Fmap_$starting_eval\.tsv w]
# # 			puts $out $fam_map_str
# # 			close $out

# # 			puts "Master counter: [expr $master_counter + 1]"
# # 			incr master_counter
# # 		}
		
# # 	}
# # 	db1 close
# # }

# # ###########################################################################

# # ###
# # set runt_output 1
# # ###
# # if {$final_groups_switch == 1} {

# # 	cd $direct/$ival
# # 	set eval_dirs [reverse_dict [glob -type d -- -*]]
# # 	set bottom_evalue [lindex $eval_dirs end]
# # 	openfile $fmap_path/Fmaps/Fmap_$bottom_evalue\.tsv
# # 	set base_families [split [string trim $data] \n]

# # 	set fmap_dirs [lrange $eval_dirs 0 end-1]
# # 	set base_fams_done {}
# # 	set runt_fams {}
# # 	set family_group_table "Final number\tE-value\tOriginal number\tSequence number\n"
# # 	file mkdir $output_path/$output_dir_name

# # 	foreach dir $fmap_dirs {
# # 		openfile $fmap_path/Fmaps/Fmap_$dir\.tsv
# # 		set fmaps [split [string trim $data] \n]

# # 		## Depending on which family map file we are looking at, we need to make sure we probe the correct evalue group fasta file. Here we set the correct range of evalue folders to look into depending on the type of fmap ##
# # 		set test_fmap_len [split [lindex $fmaps 0] \t]
# # 		set fmap_index [expr [llength $eval_dirs] - [llength $test_fmap_len]]
# # 		set group_fastas_dirs [lrange $eval_dirs $fmap_index end]
		

# # 		set i 0
# # 		foreach fmap $fmaps {
# # 			set fam_numbers [split $fmap \t]
# # 			## The families are eventually named based on the group number at the most basal (e.g. -10) level ##
# # 			set final_group_number [lindex $fam_numbers end]

# # 			## If we want to pick the best group that still contains all the Geobacillus sequences that are in the basal group, do the following... ##
# # 			if {[string tolower $minimum_geo_number] == "all"} {
# # 				set minimum_geo_seqs [regexp -all {\(Geobac\)} [openfile $direct/$ival/$bottom_evalue/Group_fastas_$ival/$final_group_number\.faa]]
# # 			} else {
# # 				set minimum_geo_seqs $minimum_geo_number
# # 			}

# # 			## Ensure we have not processed this group yet ##
# # 			if {[lsearch -exact $base_fams_done $final_group_number] == -1} {
# # 				set group_fasta_dir_counter 0

# # 				foreach group $fam_numbers {
# # 					## Make sure we are looking in the right evalue for the group number ##
# # 					set sample_evalue [lindex $group_fastas_dirs $group_fasta_dir_counter]
# # 					openfile $direct/$ival/$sample_evalue/Group_fastas_$ival/$group\.faa
# # 					set gene_no [regexp -all {>} $data]
# # 					set geobac_no [regexp -all {\(Geobac\)} $data]

# # 					# If this is not the bottom evalue stringency, and the context is greater than the desired value, and the Geobac no. is greater than the desired value, then use the family at this level, if not - continue #
# # 					if {$group_fasta_dir_counter < [expr [llength $group_fastas_dirs] - 1] && $gene_no >= $minimum_gene_context && $geobac_no >= $minimum_geo_seqs} {

# # 						file copy -force $direct/$ival/$sample_evalue/Group_fastas_$ival/$group\.faa $direct/$ival
# # 						file rename -force $direct/$ival/$group\.faa $output_path/$output_dir_name/$final_group_number\.faa

# # 						lappend base_fams_done $final_group_number
# # 						append family_group_table "$final_group_number\t$sample_evalue\t$group\t$gene_no\n"
# # 						incr i
# # 						puts "$i === $final_group_number"
# # 						break
# # 					# If this IS the bottom evalue go here #
# # 					} elseif {$group_fasta_dir_counter == [expr [llength $group_fastas_dirs] - 1]} {
# # 						# If the gene number is >= 4, a tree can be reconstructed - so we will use this family #
# # 						if {$gene_no >= 4} {

# # 							file copy -force $direct/$ival/$sample_evalue/Group_fastas_$ival/$group\.faa $direct/$ival
# # 							file rename -force $direct/$ival/$group\.faa $output_path/$output_dir_name/$final_group_number\.faa

# # 							lappend base_fams_done $final_group_number
# # 							append family_group_table "$final_group_number\t$sample_evalue\t$group\t$gene_no\n"
# # 							incr i
# # 							puts "$i === $final_group_number"
# # 							break
# # 						# If the are fewer than 4 genes in the group (even at the lowest evalue stringency e.g. -10); set the family as a runt #
# # 						} else {
# # 							lappend base_fams_done $final_group_number
# # 							lappend runt_fams "$final_group_number\t$dir"
# # 							append family_group_table "$final_group_number\t$sample_evalue\t$group\tN/A\n"
# # 							incr i
# # 							puts "$i === $final_group_number"
# # 							break
# # 						}
# # 					} else {
# # 						incr group_fasta_dir_counter
# # 					}
# # 				}
# # 			} else {
# # 				incr i
# # 				puts "$i === SKIPPED"
# # 				continue
# # 			}
# # 		}
# # 		puts "\n\nFINAL: $i / [llength $fmaps]"
# # 		incr fmap_dir_counter
# # 	}

# # 	#######

# # 	foreach group $base_families {
# # 		if {[lsearch -exact $base_fams_done $group] == -1} {
# # 			openfile $direct/$ival/$bottom_evalue/Group_fastas_$ival/$group\.faa
# # 			set gene_no [regexp -all {>} $data]
# # 			if {$gene_no >= 4} {
# # 				file copy -force $direct/$ival/$bottom_evalue/Group_fastas_$ival/$group\.faa $output_path/$output_dir_name/
# # 				lappend base_fams_done $group
# # 				append family_group_table "$group\t$bottom_evalue\t$group\t$gene_no\n"
# # 			} else {
# # 				lappend base_fams_done $group
# # 				lappend runt_fams "$group\t$bottom_evalue"
# # 				append family_group_table "$group\t$bottom_evalue\t$group\tN/A\n"
# # 			}
# # 		} else {
# # 			continue
# # 		}
# # 	}

# # 	if {$runt_output == 1} {
# # 		set out [open $output_path/Final_runt_groups.tsv w]
# # 		set runt_fams_str [join $runt_fams \n]
# # 		puts $out [string trim $runt_fams_str]
# # 		close $out
# # 	}

# # 	set out [open $output_path/Final_family_group_list.tsv w]
# # 	puts $out [string trim $family_group_table]
# # 	close $out
# # }

# # ###########################################################################

# # puts "ALL DONE!"