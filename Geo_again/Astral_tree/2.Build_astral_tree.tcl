#!/usr/local/bin/tclsh
source ~/Documents/Scripts/General_utils.tcl

set direct			/Users/aesin/Desktop/Geo_again
set gene_noDup_dir	$direct/Group_fastTree/Final_trees_noDup
set astral_out_dir	$direct/Astral_speciesTree
set astral_execute	/Users/aesin/Desktop/Software/Astral/astral.5.5.9.jar

file mkdir $astral_out_dir

## Combine input gene trees into a single master file
puts stdout	"Combining gene trees..."
set input_gene_trees	[glob $gene_noDup_dir/*txt]
set combined_tree_l		{}
set combined_tree_name	"Astral_noDup_geneTrees.txt"
foreach tree $input_gene_trees {
	set tree_data [string trim [openfile $tree]]
	lappend output_tree_l $tree_data
}
set out [open $astral_out_dir/$combined_tree_name w]
puts $out [join $output_tree_l \n]
close $out

## Run Astral
set astral_out_name	"Astral_noDup_speciesTree.txt"
chan configure stdout -buffering none
puts stdout	"Runnings Astral..."
exec >&@stdout java -Xms4g -Xmx60g -jar $astral_execute -i $astral_out_dir/$combined_tree_name -o $astral_out_dir/$astral_out_name
