## FIX TREE ##
source ~/Dropbox/Scripts/General_utils.tcl

set direct /users/aesin/desktop/Geo_v_all/2.0/Tree_build; cd $direct
set sizes [glob -type d *]

set pata {k+?[0-9]+?:+?}
set patb {\t+?.+?\n+?}

foreach size $sizes {
	cd $direct/$size
	puts "$size"
	set tree_dirs [glob -type d *]

	foreach tree_dir $tree_dirs {
		cd $direct/$size/$tree_dir
		set out_number $tree_dir

		set final_file ""
		set tree ""

		set final_file [lindex [glob -nocomplain *_tree.txt] 0]

		if {$final_file eq ""} {

			set besttrees [glob -nocomplain *bipartitions.*]
			if {[llength $besttrees] == 0} {set besttrees [glob -nocomplain *bestTree.*]}
			set tree [lindex $besttrees 0]

			set key [lindex [glob KEY*] 0]
			openfile $key; set key $data

			if {$tree ne ""} {
				openfile $tree
				set branchids [regexp -all -inline $pata $data]
				foreach k_val $branchids {
					regsub {:} $k_val {} k_val_trunc
					regexp -line $k_val_trunc$patb $key hit
					set binomial [lindex [split [string trim $hit] \t] 1]
					set class [lindex [split [string trim $hit] \t] 2]
					set id [lindex [split [string trim $hit] \t] 3]
					# If you use k_val_trunc below, you may substitute k17 when you want to find k1 - so searching for k1: is more accurate #
					regsub $k_val $data "$binomial \{$class\} \{$id\}:" data
				}

				set out [open $out_number\_tree.txt w]; puts $out $data; close $out
				continue
			} else {
				puts $tree_dir
				break
			}
		}

		set tree_data [string trim [openfile $final_file]]

		if {[regexp {@} $tree_data] != 0 || [regexp {€} $tree_data] != 0 || [regexp {\$} $tree_data] != 0} {
			puts "File $size/$final_file == corrupted... Redoing"
			file delete $final_file

			set key [lindex [glob KEY*] 0]
			openfile $key; set key $data

			set besttrees [glob -nocomplain *bipartitions.*]
			if {[llength $besttrees] == 0} {set besttrees [glob -nocomplain *bestTree.*]}
			set tree [lindex $besttrees 0]

			##########################################################################
			## Relabel the tips, and move the resultant folders/files to the output ##
			## destinations ONLY if there is a "best tree" 							##
			##########################################################################
			if {[llength $tree] > 0} {  
				openfile $tree
				set branchids [regexp -all -inline $pata $data]
				foreach k_val $branchids {
					regsub {:} $k_val {} k_val_trunc
					regexp -line $k_val_trunc$patb $key hit
					set binomial [lindex [split [string trim $hit] \t] 1]
					set class [lindex [split [string trim $hit] \t] 2]
					set id [lindex [split [string trim $hit] \t] 3]
					# If you use k_val_trunc below, you may substitute k17 when you want to find k1 - so searching for k1: is more accurate #
					regsub $k_val $data "$binomial \{$class\} \{$id\}:" data
				}

				set out [open $out_number\_tree.txt w]; puts $out $data; close $out
			}
		}
	}
}
