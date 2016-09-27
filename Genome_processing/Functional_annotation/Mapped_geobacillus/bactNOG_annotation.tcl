#!/usr/local/bin/tclsh
## Script to hmmscan all the mappable Geobacillus sequences against the bactNOG database, and identify the COG IDs ##
source ~/Dropbox/Scripts/General_utils.tcl
source ~/Dropbox/Scripts/Procs/EggNOG.tcl

###########################################################################

## Set necessary paths ##
set direct /users/aesin/desktop
set working_directory $direct/Geo_analysis/Geo_ortholog_nucl
set bactNOG_db_path $direct/EggNOG/bactNOG_db

set global_output_path $working_directory/Functional_annotation
set hmmscan_output_path $working_directory/Functional_annotation/Hmm_outputs
set prot2cog_output_path $working_directory/Functional_annotation/Prot2COG

## Make the output directories ##
file mkdir $hmmscan_output_path
file mkdir $prot2cog_output_path

## Grab the input files ##
cd $working_directory/Groups_all_map_geo_only
set group_fasta_files [lsort -dictionary [glob *faa]]

## Progress bar ##
progress_init [llength $group_fasta_files]
set counter 0

## Make the global output lists ##
set group2cog_table [list "Group\tCOG_categories\tBroad_COG"]
set global_narrow_COG_list {}
set global_broad_COG_list {}
set global_groups_no_hit {}

## Process each one in turn ##
foreach group $group_fasta_files {
	set group_number [string range $group 0 end-4]

	## Run the hmmscan ##
	if {[file exists $hmmscan_output_path/$group_number\_hmmscan_out.txt] != 1} {
		catch {exec hmmscan --cpu 6 --tblout $hmmscan_output_path/temp_out.txt --noali $bactNOG_db_path/bactNOG.txt $group}
		file rename $hmmscan_output_path/temp_out.txt $hmmscan_output_path/$group_number\_hmmscan_out.txt
	}
		
	## Process the hmmscan output and select the COG catergory for the top hit of each protein sequence ##
	set output_table [parse_nog_search_output $hmmscan_output_path/$group_number\_hmmscan_out.txt]

	## If there are no hits, parse_nog_search_output returns a list with 0 elements ##
	if {[llength $output_table] == 0} {
		lappend global_groups_no_hit $group

		set out [open $prot2cog_output_path/$group_number\_prot2cog.tsv w]
		puts $out "No hits"
		close $out

		incr counter
		progress_tick $counter

		continue
	}

	set cog_id_full_table [eggnog_cluster_to_COG_ID $output_table redundant_COG_list]

	## Count all the functional categories for a particular group and sort them in decreasing order. E.g. {O 5} {J 5} {N 2}
	set redun_COG_count [lsort -integer -index 1 -decr [lcount $redundant_COG_list]]

	## Check whether there are an equal top number of hits to a category ##
	set COG_cat_list [list [lindex [lindex $redun_COG_count 0] 0]]
	if {[llength $redun_COG_count] > 1} {
		set index 0

		## As long as the next category had an equivalent number of hits, add it to the list as well. E.g. O J ##
		while 1 {
			if {[lindex [lindex $redun_COG_count $index] 1] == [lindex [lindex $redun_COG_count [expr $index + 1]] 1]} {
				lappend COG_cat_list [lindex [lindex $redun_COG_count [expr $index + 1]] 0]
				incr index
			} else {
				break
			}
		}
	}

	## Check that any categories listed consist of a single letter, if not, split the categories and add them individually. This should not be necesary based on the output of the eggnog_cluster_to_COG_ID proc ##
	foreach category $COG_cat_list {
		if {[string length $category] > 1} {
			set multi_cats [split $category {}]
			foreach multi_cat $multi_cats {
				lappend COG_cat_list $multi_cat
			}
			set COG_cat_list [lremove $COG_cat_list $category]
			puts $COG_cat_list
		}
	}

	## Remove any redundant letters ##
	set COG_cat_list [lsort -unique $COG_cat_list]
	
	## Translate these narrow categories into a broad COG category. If one group maps to two or more narrow categories, but some of them map to the same broad category, then output the minimum number of broad cats. ##
	set COG_broad_list {}
	foreach narrow_COG $COG_cat_list {
		if {$narrow_COG eq ""} {
			puts stderr "Group: $group_number\nNarrow COGs: $COG_cat_list"
		}
		set broad_COG [dict get $broad_category_dict $narrow_COG]
		lappend COG_broad_list $broad_COG
	}
	set COG_broad_list [lsort -unique $COG_broad_list]

	## Add the values to the global tables and lists ##
	lappend group2cog_table "$group_number\t[join $COG_cat_list " | "]\t[join $COG_broad_list " | "]"
	set global_narrow_COG_list [concat $global_narrow_COG_list $COG_cat_list]
	set global_broad_COG_list [concat $global_broad_COG_list $COG_broad_list]

	## Write out the output ##
	set out [open $prot2cog_output_path/$group_number\_prot2cog.tsv w]
	puts $out [join $cog_id_full_table \n]
	close $out

	## Progress bar and stdout ##
	incr counter
	progress_tick $counter
}

## The reference table ##
set out [open $global_output_path/Group2COG_table.tsv w]
puts $out [join $group2cog_table \n]
close $out

## The list of narrow categories (for plotting) ##
set out [open $global_output_path/Narrow_COG_list.txt w]
puts $out [join $global_narrow_COG_list \n]
close $out

## The list of broad categories ##
set out [open $global_output_path/Broad_COG_list.txt w]
puts $out [join $global_broad_COG_list \n]
close $out

## The list of groups that had no hits vs database ##
set out [open $global_output_path/Groups_no_ontology.txt w]
puts $out [join $global_groups_no_hit \n]
close $out
