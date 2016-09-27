## How many trees are there in the Tree_build directory? ##

set directory /scratch/ade110/Geo_v_all/2.0

cd $directory/Tree_build
set sizes {small large giant colossal individual}

set total_done 0
set all_completed {}
foreach size $sizes {
	if {[file exists $directory/Tree_build/$size] != 1} {
		continue
	}
	cd $directory/Tree_build/$size
	set completed_directories [glob -nocomplain -type d *]
	set all_completed [concat $all_completed $completed_directories]
	set num_completed_directories [llength [glob -nocomplain -type d *]]
	set total_done [expr $total_done + $num_completed_directories]
}
puts "Total trees completed: $total_done"

# # Final set of trees to be completed #
# set to_be_completed {} 
# set n 1
# while {$n < 6050} {
# 	if {[lsearch $all_completed $n] == -1} {
# 		lappend to_be_completed $n
# 	}
# 	incr n
# }

# foreach alignment $to_be_completed {
# 	file mkdir $directory/MSA_split/Final_set/$alignment
# 	file copy $directory/Group_fastas_MSA/$alignment\.fas $directory/MSA_split/Final_set/$alignment/
# 	file copy $directory/Group_fastas_MSA/KEY_$alignment\.txt $directory/MSA_split/Final_set/$alignment/
# }

## How many trees are there left in the submission directories? ##

cd $directory/MSA_split
set core_groups [glob -type d *]
foreach core_group $core_groups {
	cd $directory/MSA_split/$core_group
	set sizes {small large giant colossal individual}

	foreach size $sizes {
		cd $directory/MSA_split/$core_group/$size
		if {$size ne "individual"} {
			
			set split_dirs [glob -type d *]
			set total_not_done 0

			foreach dir $split_dirs {
				cd $directory/MSA_split/$core_group/$size/$dir
				set not_done_trees [llength [glob -nocomplain -type d *]]
				set total_not_done [expr $total_not_done + $not_done_trees]
			}

			puts "Note done in $core_group/$size: $total_not_done"
		} else {

			set total_not_done [llength [glob -nocomplain -type d *]]
			puts "Note done in $core_group/$size: $total_not_done"

		}
	}
}

## Let's get a list of trees that are not done ##


# # OPTIONAL - move any incomplete directories for 4/8-core into the corresponding "individual" folder ##
# set groups_to_move {8_core}

# set directory /scratch/ade110/Geo_v_all/2.0

# foreach group $groups_to_move {
# 	cd $directory/MSA_split/$group
# 	set sizes {small large giant colossal}

# 	puts "$group individual trees before move: [llength [glob -nocomplain -type d $directory/MSA_split/$group/individual/*]]"

# 	set moved_counter 0

# 	foreach size $sizes {
# 		cd $directory/MSA_split/$group/$size
# 		set split_dirs [glob -type d *]

# 		foreach dir $split_dirs {
# 			cd $directory/MSA_split/$group/$size/$dir
# 			set not_done_trees [glob -nocomplain -type d *]
			
# 			foreach not_done_tree $not_done_trees {
# 				file rename -force $not_done_tree $directory/MSA_split/$group/individual
# 				incr moved_counter
# 			}
# 		}
# 	}

# 	puts "Number of tree folders moved: $moved_counter"
# 	puts "$group individual trees after move: [llength [glob -nocomplain -type d $directory/MSA_split/$group/individual/*]]"
# }

## Which individual are not done ##

cd $directory/Tree_build/individual
set ind_trees_done [glob -type d *]

cd $directory/MSA_split/4_core/individual
set small_ind_trees [glob -type d *]
cd $directory/MSA_split/8_core/individual
set large_ind_trees [glob -type d *]

set ind_trees [concat $small_ind_trees $large_ind_trees]

set ind_trees_incomplete {}
foreach ind_tree $ind_trees {
	if {[lsearch $ind_trees_done $ind_tree] == -1} {
		lappend ind_trees_incomplete $ind_tree
	}
}

puts [llength $ind_trees_incomplete]

set x 0
cd $directory/MSA_split/4_core/individual
foreach tree $small_ind_trees {
	if {[lsearch $ind_trees_incomplete $tree] == -1} {
		file delete -force $tree
	}
}

cd $directory/MSA_split/8_core/individual
foreach tree $large_ind_trees {
	if {[lsearch $ind_trees_incomplete $tree] == -1} {
		file delete -force $tree
	}
}
