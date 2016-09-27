source ~/Dropbox/Scripts/General_utils.tcl

###########################################################################

set direct /users/aesin/desktop/Geo_analysis/Geo_v_all/2.0
set tree_build_dir $direct/Tree_build_final
set nbs_dir $direct/NBS_trees
set mowgli_dir /users/aesin/desktop/Mowgli

cd $mowgli_dir/Mowgli_outputs
set mow_out_dirs [lsort -dictionary [glob -type d *out*]]

cd $tree_build_dir
set tree_dirs [glob -type d *]
puts [llength $tree_dirs]

###########################################################################
## Sort the final tree build directory (replaced the old trees with the ones rebuilt to contain all Geobacillus) into those that contain bootstraps and those that do not ##
set nbs_list {}

foreach dir $tree_dirs {
	puts $dir
	cd $tree_build_dir/$dir
	set bs_file [glob -nocomplain *bootstrap*]
	if {[llength $bs_file] == 0} {
		lappend nbs_list $dir
		file copy -force $direct/Tree_build_final/$dir $direct/Trees/Final_NBS/
	} else {
		file copy -force $dir\_tree.txt $direct/Trees/Final_BS_trees/
	}
}

puts [join [lsort -dictionary $nbs_list] \t]
puts [llength $nbs_list]


###########################################################################

## These are the trees that have been rebuilt to include all Geobacillus ##
cd $direct/Tree_build_no_missing_geo
set all_geo_dirs [lsort -dictionary [glob -type d *]]

foreach mow_out_dir $mow_out_dirs {
	cd $mowgli_dir/Mowgli_outputs/$mow_out_dir
	file mkdir $mowgli_dir/OLD/Mowgli_outputs/$mow_out_dir
	set present 0

	foreach all_geo_dir $all_geo_dirs {
		if {[file exists $mowgli_dir/Mowgli_outputs/$mow_out_dir/$all_geo_dir] == 1} {
			incr present
			file rename $mowgli_dir/Mowgli_outputs/$mow_out_dir/$all_geo_dir $mow_out_dir/OLD/Mowgli_outputs/$mow_out_dir/$all_geo_dir
		}
	}

	puts "For $mow_out_dir there are $present folders to be deleted"
}



###########################################################################

## Find how many still need to be reconciled ##

set one_mow_out_dir [lindex $mow_out_dirs 0]
cd $mowgli_dir/Mowgli_outputs/$one_mow_out_dir
set reconciled_done [glob -type d *]

cd $direct/Tree_build_final
set all_built_trees [glob -type d *]

set not_reconciled {}
foreach built_tree $all_built_trees {
	puts $built_tree
	if {[lsearch $reconciled_done $built_tree] == -1} {
		lappend not_reconciled $built_tree
	}
}

foreach not_recon $not_reconciled {
	if {[file exists $direct/Trees/Final_NBS/$not_recon] == 0} {
		file copy $direct/Trees/Final_BS_trees/$not_recon\_tree.txt $mowgli_dir/Reconcile_BS/
	}
}
