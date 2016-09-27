#!/usr/local/bin/tclsh

### Get functional annotation for Anoxy/Geo and Geo only proteins ###

source ~/Dropbox/Scripts/General_utils.tcl
package require math::statistics

###########################################################################

set direct /users/aesin/desktop/Geo_analysis/Geo_v_all/2.0/Geo_AG_only_groups
file mkdir $direct/Combined_map_fastas

set not_all_geo_map [split [string trim [openfile /users/aesin/desktop/Geo_analysis/Geo_ortholog_nucl/Removed_groups.tsv]] \n]
set group2cog_tab [split [string trim [openfile /users/aesin/desktop/Geo_analysis/Geo_ortholog_nucl/Functional_annotation/Group2COG_table.tsv]] \n]

###########################################################################
set x 0
set do_not_map 0
set poor_char 0
set have_COG {}

## Get all the group numbers ##
set group_numbers {}
cd $direct/Combined_fastas
set fastas [lsort -dictionary [glob *faa]]
foreach file $fastas {
	set g_number [string range $file 0 end-4]

	if {[lsearch $not_all_geo_map $g_number] != -1} {
		incr do_not_map
	} else {
		file copy -force $file $direct/Combined_map_fastas
	}

	set line [lsearch -glob -inline $group2cog_tab $g_number*]
	if {[string length $line] != 0} {
		incr x
		if {[lindex [split $line \t] 2] eq "Poorly Characterized"} {
			incr poor_char
		} else {
			lappend have_COG $line
		}
	}

	
}

puts "$x\tPC: $poor_char"
puts $do_not_map
puts "\n[join $have_COG \n]"
