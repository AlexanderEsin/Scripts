#!/usr/local/bin/tclsh
set folder All_fastas_ncbi_id
set file_ext "_protein.faa"
set direct /users/aesin/desktop/Clean_Proteomes
set output_table_switch 1


## Source procs ##
source ~/Dropbox/Scripts/General_utils.tcl

proc taxa_name_clean {name} {
	global newname
	regsub -all "/" $name "_" newname
	regsub -all -- {\-} $newname "_" newname
	regsub -all {=} $newname {} newname
	regsub -all "  " $newname " " newname
	regsub -all " " $newname "_" newname
	return $newname
}

###########################################################################

## Open the assembly reference file. It is IMPORTANT to use the old one to pick up all the NCBI IDs correctly##
cd ~/desktop/Assembly_summaries
openfile Oct14_assembly_summary.txt
set reflist $data

cd $direct/$folder
set orgs [glob *$file_ext]

set translation_table {}
set pata {.+?_+?.+?\.+?[0-9]+?}
set patb {.+?\n+?}
set i 1
set num_missed 0

foreach org $orgs {
	regsub "$file_ext" $org {} org
	regexp -line $pata $org hit

	set hit_in_reflist [regexp -line $hit$patb $reflist match]
	if {$hit_in_reflist == 0} {
		puts "No match in reflist for: $org"
		incr num_missed
	}

	set match [string trim $match]
	set match [split $match \t]
	set newname [lindex $match 7]

	
	set newname [taxa_name_clean $newname]
	# file rename -force $org$file_ext $newname$file_ext
	lappend translation_table "$org\t$newname"

	puts "Translated ... $i/[llength $orgs]"
	incr i
}

if {$output_table_switch == 1} {
	set out [open $direct/$folder\_binomial_translation.tsv w]
	puts $out [join $translation_table \n]
	close $out
}

puts "Not Translated: $num_missed"

######