#!/usr/local/bin/tclsh
## Funtional annotation of consistent Mowgli results true for horizontal transfer ##

source ~/Dropbox/Scripts/General_utils.tcl
source ~/Dropbox/Scripts/Procs/EggNOG.tcl

###########################################################################
set inside_included TRUE
set all_scenarios FALSE

set direct /users/aesin/desktop
set input_direct $direct/Mowgli/Consistent_HGT_Vertical
set output_direct $direct/Geo_analysis/HGT_functional_annotation/Consistent_HGT

###############################

if {$inside_included == TRUE} {
	if {$all_scenarios == TRUE} {
		set input_path  $input_direct/Consistent_Full/All_scenarios
		set output_direct $output_direct/Full/All_scenarios
	} else {
		set input_path  $input_direct/Consistent_Full/Scenarios_1_2
		set output_direct $output_direct/Full/Scenarios_1_2
	}
} else {
	if {$all_scenarios == TRUE} {
		set input_path  $input_direct/Consistent_NO_IG/All_scenarios
		set output_direct $output_direct/NO_IG/All_scenarios
	} else {
		set input_path  $input_direct/Consistent_NO_IG/Scenarios_1_2
		set output_direct $output_direct/NO_IG/Scenarios_1_2
	}
}

###########################################################################

## Open the functional annotation table ##
set func_annot_table [split [string trim [openfile $direct/Geo_analysis/Geo_ortholog_nucl/Functional_annotation/Group2COG_table.tsv]] \n]

## Removed groups list ##
set removed_list [split [string trim [openfile $direct/Geo_analysis/Geo_ortholog_nucl/Removed_groups.tsv]] \n]

cd $input_path
set input_files [glob true_*]

###############################

foreach consistent_result $input_files {

	## Is this result horizontal or vertical? ##
	if {[string match "*HGT*" $consistent_result] == 1} {
		set type HGT
	} else {
		set type Vertical
	}

	## Penalties being tested ##
	regsub -all "t" [string range [regexp -inline {t+[0-9]+.+(?:.tsv)+} $consistent_result] 0 end-4] {} penalty_str

	#puts "$consistent_result\t$type"

	set consistent_data [split [string trim [openfile $consistent_result]] \n]

	###############################

	set functional_list {}
	set broad_COG_list {}
	set groups_not_found_l {}
	set translation_hgt {}

	set groups_counter 0
	set no_func_counter 0

	###############################

	foreach group $consistent_data {
		
		set functional_annot_line ""
		set functional_annot_line [lsearch -glob -inline $func_annot_table $group\t*]

		if {$functional_annot_line eq ""} {
			# puts "Fucntional annotation not done: $group"
			incr no_func_counter
			lappend groups_not_found_l $group
			
			continue
		}

		set temp_broad_COG_list {}

		## Some groups may have more than one COG assigned. We want to add all of them in ##
		set functional_group [lindex [split $functional_annot_line \t] 1]
		set functional_groups [split $functional_group \|]

		## For each potential COG.. ##
		foreach annot $functional_groups {
			set annot [string trim $annot]

			## Add each COG to the list of COGs ##
			lappend functional_list $annot
			## If it's a J and is from the HGT set, add it to a list ##
			if {$type eq "HGT" && $annot eq "J"} {
				lappend translation_hgt $group
			}
			
			## If two or more COGs are assigned, they can have different "broad" categories ##
			set broad_COG [dict get $broad_category_dict $annot]
			lappend temp_broad_COG_list $broad_COG
		}

		## Add any -different- broad categories to the global list ##
		set temp_broad_COG_list [lsort -unique $temp_broad_COG_list]
		foreach broad_cog $temp_broad_COG_list {
			lappend broad_COG_list $broad_cog
		}

		incr groups_counter
	}

	set output_path $output_direct/$type/$penalty_str
	file mkdir $output_path

	set narrow_out [open $output_path/Narrow_COG.tsv w]
	puts $narrow_out [join $functional_list \n]
	close $narrow_out

	set broad_out [open $output_path/Broad_COG.tsv w]
	puts $broad_out [join $broad_COG_list \n]
	close $broad_out

	if {$type eq "HGT"} {
		set cogJ_out [open $output_path/J_HGT_list.txt w]
		puts $cogJ_out [join $translation_hgt \n]
		close $cogJ_out
	}
	

	set log_out [open $output_path/Stats.txt w]
	puts $log_out "Consistent set: $type\nPenalties tested: $penalty_str\nGroups tested: $groups_counter\nNo functional annotation found: $no_func_counter"
	close $log_out

	puts stdout "\n###############################\nConsistent set: $type\nPenalties tested: $penalty_str\nGroups tested: $groups_counter\nNo functional annotation found: $no_func_counter"
}



