#!/usr/local/bin/tclsh
source ~/Dropbox/Scripts/General_utils.tcl
package require sqlite3

## Set paths ##
set direct /users/aesin/desktop
set working_directory $direct/Geo_analysis/Geo_ortholog_nucl
set mowgli_results_directory $direct/Mowgli/Mowgli_outputs

set removed_groups_list [split [string trim [openfile $working_directory/Removed_groups.tsv]] \n]
set mowgli_results [split [string trim [openfile $mowgli_results_directory/T3_D2_L1_results.txt]] \n]

set reconciled 0
set not_reconciled 0
set vertical 0
set horizontal 0
set to_be_reconc_list {}

foreach removed_group $removed_groups_list {
	set hgt_into_anoxy_geo [lindex [split [lsearch -inline -glob $mowgli_results $removed_group\t*] \t] 1]
	if {[string length $hgt_into_anoxy_geo] != 0} {
		incr reconciled
		if {$hgt_into_anoxy_geo == 0} {
			incr vertical
		} elseif {$hgt_into_anoxy_geo == 1} {
			incr horizontal
		} else {
			puts stderr "Should not have a non-binary result"
			exit 2
		}
	} else {
		incr not_reconciled
		lappend to_be_reconc_list $removed_group
	}
}




set out [open $working_directory/Removed_groups_vertical_horizontal.tsv w]
puts $out "Set\tType\tNumber"
puts $out "Removed\tVertical\t$vertical"
puts $out "Removed\tHorizontal\t$horizontal"
## These numbers are based on the IG-included t3 reconciliation ##
puts $out "Total\tVertical\t1921"
puts $out "Total\tHorizontal\t3647"
close $out

# puts "Number of removed groups: [llength $removed_groups_list]"
# puts "Number of those that were reconciled by Mowgli (t3): $reconciled"
# puts "\tOf those, vertical signal:	$vertical"
# puts "\tand horizontal signal:		$horizontal"

puts [join $to_be_reconc_list \n]