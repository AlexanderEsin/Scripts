#!/usr/local/bin/tclsh
source ~/Documents/Scripts/General_utils.tcl

set direct		"/Users/aesin/Desktop/Geo_again/Functional_annotation/Output_annots/"
set annot_files	[glob $direct/*emapper.annotations]

foreach file $annot_files {
	set data	[string trim [openfile $file]]
	regsub -all {#query} $data {query} new_data
	set out		[open $file w]
	puts $out	$new_data
	close $out
}