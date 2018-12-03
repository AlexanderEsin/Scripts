#!/usr/bin/tclsh

source /home/ade110/Scripts/General_utils.tcl

set eval_cutoff	[lindex $argv 0]
set eval_num	[string range $eval_cutoff 2 end]

## Directories and file paths
set direct			/home/ade110/Work/Bacillus/RBBH
set master_dir		$direct/Master_RBBH
set master_file		$master_dir/Master_RBBH_$eval_num\.txt
set new_master		$master_dir/Master_RBBH_adj_$eval_num\.txt

## Exit if the original master file does not exist
if {[file exists $master_file] != 1} {
	puts stderr "Cannot find file $master_file"
	exit
}

## Open the adjusted master file for appending
set newOut		[open $new_master a]
## Open old master file for reading
set fparse		[open $master_file r]

## Line by line change any 0.0 evalues to 1e-200
set line_count	0
while {[gets $fparse data] >= 0} {
	incr line_count
	set line_split		[split $data \t]
	if {[lindex $line_split end] eq "0.0"} {
		set new_line	[join [concat [lrange $line_split 0 1] [list "1e-200"]] \t]
		puts $newOut	$new_line 
	} else {
		puts $newOut	$data
	}
}

## Close the file channels
close $fparse
close $newOut

## Check the number of lines in the new master file
set new_line_count	0
set newOut	[open $new_master r]
while {[gets $newOut data] >= 0} {
	incr new_line_count
}

puts stdout	"Old master file at eval : $eval_num had $line_count lines"
puts stdout	"New master file at eval : $eval_num has $new_line_count lines"