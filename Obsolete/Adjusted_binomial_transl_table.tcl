#!/usr/local/bin/tclsh

# The Prok_gbff_files_ID_binomial_translation.tsv used to have ncbi ids suffixed with "_genomic". Wanted to remove it as it is not necessary. Complementary changes made to both the "Rename_gbff_binomial.tcl" and "Calc_gc_content_fam.tcl" scripts #


## Source procs ##
source ~/Dropbox/Scripts/General_utils.tcl

set direct /users/aesin/desktop/Clean_Genomes

cd $direct
set table [string trim [openfile Prok_gbff_files_ID_binomial_translation.tsv]]
set entries [split $table \n]
puts "Number of entries: [llength $entries]"
set new_table {}

foreach entry $entries {
	set old_id [lindex [split $entry \t] 0]
	regsub {_genomic} $old_id {} adj_id
	set new_entry "$adj_id\t[lindex [split $entry \t] 1]"
	lappend new_table $new_entry
}

puts "Number of adjusted entries: [llength $new_table]"
set new_table [join $new_table \n]
set out [open $direct/Prok_gbff_files_ID_binomial_translation_temp.tsv w]
puts $out $new_table
close $out