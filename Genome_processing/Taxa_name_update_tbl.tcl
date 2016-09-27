#!/usr/local/bin/tclsh
source ~/Dropbox/Scripts/General_utils.tcl
source ~/Dropbox/Scripts/Procs/gc_content.tcl

proc taxa_name_clean {name} {
	global newname
	regsub -all "/" $name "_" newname
	regsub -all -- {\-} $newname "_" newname
	regsub -all {=} $newname {} newname
	regsub -all "  " $newname " " newname
	regsub -all " " $newname "_" newname
	return $newname
}

## Combine the taxon name update since release 60 (the one I used for the Geo_v_all dataset is around 66/67). The original files are stored in Clean_Genomes/Taxon_name_update.

set direct /users/aesin/desktop/Clean_Genomes/Taxon_name_update
cd $direct
puts [exec pwd]
set update_files [glob *update]
set i 0

set taxa_name_update_new_t {}

foreach file $update_files {
	set data [string trim [openfile $file]]
	## Trim down to just the nomenclature changes ##
	regsub -all "\n\n" $data {£} data
	set changes_list [split [string trim [string range $data [string last £ $data]+1 end]] \n]
	foreach change $changes_list {
		set entries [split $change |]
		set new_name [taxa_name_clean [lindex $entries 1]]
		set old_name [taxa_name_clean [lindex $entries 3]]
		lappend taxa_name_update_new_t "$old_name\t$new_name"
	}
	if {$i == 0} {
		puts $taxa_name_update_new_t
	}
	incr i
}

set out [open /users/aesin/desktop/Clean_Genomes/Old_to_new_taxa_name_tbl.tsv w]
puts $out [join $taxa_name_update_new_t \n]
close $out