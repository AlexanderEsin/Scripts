source ~/Dropbox/Scripts/General_utils.tcl
## Combine all the newick trees in one directory, ending with extension, into a single file ##
proc combine_tree_files {directory extension output_name} {
	cd $directory
	set tree_files [glob *$extension]
	puts "Number of tree files: [llength $tree_files]"

	set output_tree_l {}
	foreach tree $tree_files {
		set tree_data [string trim [openfile $tree]]
		lappend output_tree_l $tree_data
	}

	set out [open $output_name w]
	puts $out [join $output_tree_l \n]
	close $out
}