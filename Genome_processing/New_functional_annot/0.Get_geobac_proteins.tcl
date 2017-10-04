#!/usr/local/bin/tclsh
source ~/Documents/Scripts/General_utils.tcl

# The group fasta files were first compared to the set used for trees for consistency - no difference found
# diff -r /Users/aesin/Desktop/Geo_analysis/Geo_ortholog_nucl/Groups /Users/aesin/Desktop/Geo_analysis/Geo_v_all/2.0/Final_family_groups

set master_dir		/Users/aesin/Desktop/Geo_analysis/Geo_ortholog_nucl
set full_group_dir	Groups
set func_annot_dir	New_functional_annot

set full_groups_files	[glob $master_dir/$full_group_dir/*faa]
puts "Processing [llength $full_groups_files] group files...\n"

# For each group fasta file, extract just the Geobacillus sequence. We will run the COG detection on the Geobacillus proteins only
set i 1
foreach group_file $full_groups_files {

	set file_name [string trim [file tail $group_file]]
	set all_prots	[split_genes [string trim [openfile $group_file]]]

	set geobac_prots [lsearch -all -inline $all_prots "* \(Geobac\) *"]
	set out [open $master_dir/$func_annot_dir/Input_groups/$file_name w]
	puts $out [join $geobac_prots ""]
	close $out

	puts "Done: $i / [llength $full_groups_files]"
	incr i
}
