#!/usr/local/bin/tclsh
## Functional annotation of Mowgli HGT results per penalty ##

source ~/Dropbox/Scripts/General_utils.tcl
source ~/Dropbox/Scripts/Procs/EggNOG.tcl

###########################################################################
set inside_included TRUE
set all_scenarios FALSE

set direct /users/aesin/desktop
set input_direct	$direct/Mowgli/Mowgli_outputs
set output_direct	$direct/Geo_analysis/HGT_Functional_annotation

if {$inside_included == TRUE} {
	set input_path  $input_direct/Results_Full
	if {$all_scenarios == TRUE} {
		set output_direct $output_direct/Per_penalty/Full/All_scenarios
	} else {
		set output_direct $output_direct/Per_penalty/Full/Scenarios_1_2
	}
} else {
	set input_path  $input_direct/Results_NO_IG
	if {$all_scenarios == TRUE} {
		set output_direct $output_direct/Per_penalty/NO_IG/All_scenarios
	} else {
		set output_direct $output_direct/Per_penalty/NO_IG/Scenarios_1_2
	}
}


###########################################################################
## Open the functional annotation table ##
set func_annot_table [split [string trim [openfile $direct/Geo_analysis/Geo_ortholog_nucl/Functional_annotation/Group2COG_table.tsv]] \n]

## Removed groups list ##
set removed_list [split [string trim [openfile $direct/Geo_analysis/Geo_ortholog_nucl/Removed_groups.tsv]] \n]

cd $input_path
set results_files [lsort -dictionary [glob *results.txt]]

###############################

foreach result $results_files {
	## Get the transfer penalty value ##
	set penalty [string range $result 1 [string first \_ $result]-1]

	## Make the output path ##
	set output_path "$output_direct/T$penalty"
	file mkdir $output_path

	## Get the Mowgli results ##
	set results_data [split [string trim [openfile $result]] \n]

	set hor_functional_l {}
	set ver_functional_l {}

	set hor_broad_COG_l {}
	set ver_broad_COG_l {}

	set groups_not_found_l {}
	set translation_hgt	{}

	set total 0
	set ver_counter 0
	set hor_counter 0
	set hor_no_func 0
	set ver_no_func 0
	set scenario_3_4 0

	###############################

	foreach reconcilitation $results_data {

		set group [lindex [split $reconcilitation \t] 0]
		## Type here can be 0 (vertical) or 1 (horizontal) ##
		set type [lindex [split $reconcilitation \t] 1]
		## 1, 2, 3, or 4 ##
		set scenario [lindex [split $reconcilitation \t] 2]

		###############################
		## If the we are only taking scenarios 1 & 2, ignore any other results ##
		if {$all_scenarios != TRUE} {
			if {$scenario > 2} {
				incr scenario_3_4
				incr total
				continue
			}
		}

		###############################

		set functional_annot_line ""
		set functional_annot_line [lsearch -glob -inline $func_annot_table $group\t*]

		###############################
		## If there's no functional annotation for that group... ##
		if {$functional_annot_line eq ""} {
			#puts "Fucntional annotation not done: $group"
			if {$type == 1} {
				incr hor_counter
				incr hor_no_func
			} else {
				incr ver_counter
				incr ver_no_func
			}

			lappend groups_not_found_l $group
			
			incr total
			continue
		}

		###############################
		set temp_broad_COG_list {}

		## Some groups may have more than one COG assigned. We want to add all of them in ##
		set functional_group [lindex [split $functional_annot_line \t] 1]
		set functional_groups [split $functional_group \|]

		## For each potential COG.. ##
		foreach annot $functional_groups {
			set annot [string trim $annot]

			if {$type == 1} {
				## Add each COG to the list of COGs ##
				lappend hor_functional_l $annot
				## If it's a J and is from the HGT set, add it to a list ##
				if {$annot eq "J"} {
					lappend translation_hgt $group
				}
			} else {
				lappend ver_functional_l $annot
			}
			
			## If two or more COGs are assigned, they can have different "broad" categories ##
			set broad_COG [dict get $broad_category_dict $annot]
			lappend temp_broad_COG_list $broad_COG
		}

		## Add any -different- broad categories to the global lists ##
		set temp_broad_COG_list [lsort -unique $temp_broad_COG_list]
		foreach broad_cog $temp_broad_COG_list {
			if {$type == 1} {
				lappend hor_broad_COG_l $broad_cog
			} else {
				lappend ver_broad_COG_l $broad_cog
			}
		}

		## Increment relevant counters ##
		if {$type == 1} {
			incr hor_counter
		} else {
			incr ver_counter
		}

		incr total
	}

	###########################################################################

	## Output the results ##
	set hgt_narrow_cog_out [open $output_path/HGT_narrow_COG.tsv w]
	puts $hgt_narrow_cog_out [join $hor_functional_l \n]
	close $hgt_narrow_cog_out

	set hgt_broad_cog_out [open $output_path/HGT_broad_COG.tsv w]
	puts $hgt_broad_cog_out [join $hor_broad_COG_l \n]
	close $hgt_broad_cog_out

	set vert_narrow_cog_out [open $output_path/Vertical_narrow_COG.tsv w]
	puts $vert_narrow_cog_out [join $ver_functional_l \n]
	close $vert_narrow_cog_out

	set vert_broad_cog_out [open $output_path/Vertical_broad_COG.tsv w]
	puts $vert_broad_cog_out [join $ver_broad_COG_l \n]
	close $vert_broad_cog_out

	set hgt_cogJ_out [open $output_path/Translation_J_HGT_list.txt w]
	puts $hgt_cogJ_out [join $translation_hgt \n]
	close $hgt_cogJ_out

	# Output the stats file and the same to stdout ##
	puts stdout "\n###############################\nPenalty tested: $penalty\nTotal: $total\nVertical: $ver_counter\nHGT found: $hor_counter\nHGT not found: $hor_no_func\nVertical not found: $ver_no_func"

	set log_out [open $output_path/Stats.txt w]
	puts $log_out "Penalty tested: $penalty\nTotal: $total\nVertical: $ver_counter\nHGT found: $hor_counter\nHGT not found: $hor_no_func\nVertical not found: $ver_no_func"	
	if {$all_scenarios != TRUE} {
		puts $log_out "Non 1 or 2 scenarios removed: $scenario_3_4"
		puts stdout "Non 1 or 2 scenarios removed: $scenario_3_4"
	}
	close $log_out

}
