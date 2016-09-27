cd ~/desktop/Test_restart/1286

proc tcl::mathfunc::roundto {value decimalplaces} {
	expr {round(10**$decimalplaces*$value)/10.0**$decimalplaces}
}


## Find any incompleted Raxml intermediate ML tree files ##
set ML_result_trees [glob -nocomplain RAxML_result*]
## As long as there are at least two files, we can figure out how long it takes to build one ML tree. We cannot use the last result file as it is continuously updated and its timestamp may not represent the completed tree. We use the penultimate completed results tree - hence we need at least two ##
if {[llength $ML_result_trees] > 1} {
	set number_ML_completed [expr [llength [split [exec ls -lt | grep result] \n]] -1]
	set first_parsimony_file [lindex [split [exec ls -lt | grep parsimony | grep RUN.0]] end]; set first_parsimony_file_t [exec stat -c%Y $first_parsimony_file]
	set last_tree_file [lindex [split [lindex [split [exec ls -lt | grep result] \n] 1] " "] end]; set last_tree_file_t [exec stat -c%Y $last_tree_file]
	set total_time_taken [expr ([exec stat -c%Y $last_tree_file] - [exec stat -c%Y $first_parsimony_file]) / double(3600)]
	set time_per_tree [expr roundto([expr $total_time_taken / $number_ML_completed], 2)]
}


## For large or giant trees (or rerun of a small tree which already had a bootstrap file) we wanted 20 distinct ML trees that could be used to find the best one ##
if {$size == "small" || $size == "large" || $size == "giant"} {
	## If building 20 trees would exceed the soft walltime limit (- 3 hours for variation in run time and final optimisation), we want to reduce the number of trees to be built ##
	if {[expr $time_per_tree * 20] > [expr $walltime_limit - 3]} {
		## If we can build 15 trees in the super-soft walltime, set the number to 15 ##
		if {[expr $time_per_tree * 15] < [expr $walltime_limit - 3]} {
			set ML_tree_build_no 15
		## If we can build 10 trees in the super-soft walltime, set the number to 10 ##
		} elseif {[expr $time_per_tree * 10] < [expr $walltime_limit - 3]} {
			set ML_tree_build_no 10
		## If we cannot build at least 10 trees - we shift this folder into the individual folder for each group ##
		} else {
			cd $direct/MSA_split/$group/$size/[lindex $argv 0]
			file rename -force $directory $direct/MSA_split/$group/individual
			puts "This alignment: $group/$size/$directory takes $time_per_tree hours to construct a single ML tree and the full analysis not fit into the walltime limit: $walltime_limit... Moved to $group/individual\n"
			continue
		}
	}
## For colossal trees we wanted at least 10 distinct ML trees ##
} elseif {$size == "colossal"} {
	if {[expr $time_per_tree * 10] > [expr $walltime_limit - 3]} {
		cd $direct/MSA_split/$group/$size/[lindex $argv 0]
		file rename -force $directory $direct/MSA_split/$group/individual
		puts "This alignment: $group/$size/$directory takes $time_per_tree hours to construct a single ML tree and the full analysis not fit into the walltime limit: $walltime_limit... Moved to $group/individual\n"
		continue
}
