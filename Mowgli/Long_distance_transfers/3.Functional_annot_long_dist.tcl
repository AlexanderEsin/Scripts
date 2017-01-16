#!/usr/local/bin/tclsh
## Get the functional annotation for the refined groups list ##

source ~/Dropbox/Scripts/General_utils.tcl
source ~/Dropbox/Scripts/Procs/EggNOG.tcl

###########################################################################
## List of penalties for which to derive the COGs ##
set test_penalty_list [list 4 5 6]

## Other options ##
set inside_included TRUE
set all_scenarios FALSE

set direct /users/aesin/desktop
set input_refined_direct	$direct/Mowgli/Long_distance_HGT/
set output_func_direct		$direct/Geo_analysis/HGT_Functional_annotation/Long_distance_HGT/

if {$inside_included == TRUE} {
	if {$all_scenarios == TRUE} {
		set path  Full/All_scenarios
	} else {
		set path  Full/Scenarios_1_2
	}
} else {
	if {$all_scenarios == TRUE} {
		set path  NO_IG/All_scenarios
	} else {
		set path  NO_IG/Scenarios_1_2
	}
}

set input_refined_direct $input_refined_direct$path/Per_penalty_refined
set output_func_direct $output_func_direct$path

file mkdir $output_func_direct/Narrow; file mkdir $output_func_direct/Broad; file mkdir $output_func_direct/J_HGT; file mkdir $output_func_direct/Stats

###############################
## Open the functional annotation table ##
set func_annot_table [split [string trim [openfile /users/aesin/desktop/Geo_analysis/Geo_ortholog_nucl/Functional_annotation/Group2COG_table.tsv]] \n]

###############################
foreach penalty $test_penalty_list {

	## Get the refined results ##
	set list_of_groups [split [string trim [openfile $input_refined_direct/T$penalty\_refined_groups.tsv]] \n]
	#puts [llength $list_of_groups]

	set groups_not_found_l {}
	set functional_annot_list {}
	set broad_cog_list {}
	set translation_hgt {}

	set groups_counter 0
	set no_func_counter 0

	foreach group $list_of_groups {
		set functional_annot_line ""
		set functional_annot_line [lsearch -glob -inline $func_annot_table $group\t*]

		if {$functional_annot_line eq ""} {
			#puts "Fucntional annotation not done: $group"
			lappend groups_not_found_l $group
			
			incr groups_counter
			incr no_func_counter
			continue
		}

		set temp_broad_COG_list {}

		set functional_group [lindex [split $functional_annot_line \t] 1]
		set functional_groups [split $functional_group \|]
		foreach annot $functional_groups {
			set annot [string trim $annot]

			if {$annot eq "J"} {
				lappend translation_hgt $group
			}

			lappend functional_annot_list $annot
			
			set broad_COG [dict get $broad_category_dict $annot]
			lappend temp_broad_COG_list $broad_COG
		}

		set temp_broad_COG_list [lsort -unique $temp_broad_COG_list]
		foreach broad_cog $temp_broad_COG_list {
			lappend broad_cog_list $broad_cog
		}

		incr groups_counter
	}

	set out [open $output_func_direct/Narrow/T$penalty\_narrow_COG.tsv w]
	puts $out [join $functional_annot_list \n]
	close $out

	set out [open $output_func_direct/Broad/T$penalty\_broad_COG.tsv w]
	puts $out [join $broad_cog_list \n]
	close $out

	set hgt_cogJ_out [open $output_func_direct/J_HGT/T$penalty\_J_HGT_list.txt w]
	puts $hgt_cogJ_out [join $translation_hgt \n]
	close $hgt_cogJ_out

	set log_out [open $output_func_direct/Stats/T$penalty\_Stats.txt w]
	puts $log_out "Penalty tested: $penalty\nGroups tested: $groups_counter\nNo functional annotation found: $no_func_counter\nNumber of groups with annotation: [expr $groups_counter - $no_func_counter]"
	close $log_out

	puts stdout "\n###############################\nPenalty tested: $penalty\nGroups tested: $groups_counter\nNo functional annotation found: $no_func_counter\nNumber of groups with annotation: [expr $groups_counter - $no_func_counter]"
}

