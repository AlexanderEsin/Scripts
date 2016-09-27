###### CORRECT THE FINAL TREES (i.e. exchanging of k_values for species labels) because of the k1 vs k1: substitution problem ######

###
set direct /users/aesin/desktop/Geo_v_all
###

proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}
proc reverse_dict {lst} {
	global revdict
	while {[llength $lst] > 0} {
		set element [lindex $lst end]
		lappend revdict $element
		set lst [lrange $lst 0 end-1]
	}
}
###########################################################################

cd $direct/Geo_arch_analysis/Tree_build

set directories [glob -nocomplain -type d *]

foreach directory $directories {
	set key "KEY_$directory\.txt"

	cd $direct/Geo_arch_analysis/Geobac_archaea_aligns/outgroup
	set filelist [glob *.txt]
	if {[regexp $key $filelist] == 1} {
		openfile $key
		set key $data
		set out_flag outgroup
	} else {
		cd $direct/Geo_arch_analysis/Geobac_archaea_aligns/no_outgroup
		openfile $key
		set key $data
		set out_flag no_outgroup
	}

	cd $direct/Geo_arch_analysis/Tree_build/$directory

	set pata {k+?[0-9]+?:+?}
	set patb {\t+?.+?\n+?}

	set besttrees [glob -nocomplain *bipartitions.*]
	if {[llength $besttrees] > 0} {
		foreach tree $besttrees {
			openfile $tree
			set tree_data $data
			set branchids [regexp -all -inline $pata $tree_data]
			foreach k_val $branchids {
				regsub {:} $k_val {} k_val_trunc
				regexp -line $k_val_trunc$patb $key hit
				set binomial [lindex [split [string trim $hit] \t] 1]
				set id [lindex [split [string trim $hit] \t] 2]
				regsub $k_val $tree_data "$binomial \{$id\}:" tree_data
			}
			set out [open $direct/Geo_arch_analysis/Final_trees/$out_flag/$directory\_tree.txt w]
			puts $out $tree_data
			close $out
		}
	}
}
