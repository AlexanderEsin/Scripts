#!/usr/local/bin/tclsh
source ~/Documents/Scripts/General_utils.tcl

set files [glob "/Users/aesin/Desktop/Geo_analysis/Geo_ortholog_nucl/Functional_annot/Output_annots/*emapper.annotations"]

foreach file $files {
	set data [openfile $file]
	regsub -all {#query} $data {query} new_data
	set out [open $file w]
	puts $out $new_data
	close $out
}