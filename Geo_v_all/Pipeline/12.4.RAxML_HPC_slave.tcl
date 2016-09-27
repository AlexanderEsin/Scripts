## The ival at which the clusteres were made ##
set ival 2.0
## A soft (-1 hour) walltime limit for each run ##
set walltime_limit 71
#set direct /users/aesin/desktop/Geo_v_all/$ival
set direct /scratch/ade110/Geo_v_all/$ival
set cores 4
set group 4_core
set sizes {small large giant colossal}

##################
## Source procs ##
##################

source ~/Scripts/General_utils.tcl

###########################################################################

#######################################
## Make the necessary output folders ##
#######################################
cd $direct
if {[llength [glob -nocomplain -type d $direct/Trees/*]] == 0} {
	foreach size $sizes {
		file mkdir Trees/$size
		file mkdir Tree_build/$size
	}
}

#########################################
## Patterns for leaf name substitution ##
#########################################
set pata {k+?[0-9]+?:+?}
set patb {\t+?.+?\n+?}

########################
## Set run start time ##
########################
set start_time [clock seconds]
set soft_end_time [expr round($start_time + ($walltime_limit * 60 * 60))]

######################################
## Set the runtime start timer here ##
######################################
set current_run_start 0; set current_run_finish 0

####################################################################
## Foreach size folder, run through all the alignment directories ##
####################################################################
foreach size $sizes {

	####################################################################
	## Go into each split directory - given as argument to tcl script ##
	## and within that get all the incomplete group alignment folders ##
	####################################################################
	cd $direct/MSA_split/$group/$size/[lindex $argv 0]
	puts "Current folder directory: [exec pwd]"

	set dirs [glob -nocomplain -type d *]
	if {[llength $dirs] == 0} {
		continue
	}
	set dirs [reverse_dict $dirs]
	puts "$size === list of directories: $dirs"

	##################################
	## For each alignment directory ##
	##################################
	foreach directory $dirs {
		puts "Working on $directory ...."

		######################################################################
		## If there is already a folder in Tree_build for this directory 	##
		## and that folder has a completed final tree file, then continue	##
		## and remove the folder from the MSA_split. If the directory does 	##
		## exist in Tree_build, but there is no final file with relabelled 	##
		## tips, then DELETE the directory in Tree_build and repeat the  	##
		## reconstruction.													##
		######################################################################
		set directory_already_done_switch 0
		if {[file exists $direct/Tree_build/$size/$directory/$directory\_tree.txt] == 1} {
			puts "$size === $directory is already DONE --> removed folder from MSA"
			set directory_already_done_switch 1
			file delete -force $direct/MSA_split/$group/$size/[lindex $argv 0]/$directory
		} elseif {[file exists $direct/Tree_build/$size/$directory] == 1 && [file exists $direct/Tree_build/$size/$directory/$directory\_tree.txt] == -1} {
			puts "$size === $directory exists in Tree_build, but cannot find final file --> deleting the incomplete directory from Tree_build"
			file delete -force $direct/Tree_build/$size/$directory
		}

		if {[file exists $direct/MSA_split/$group/individual/$directory] == 1} {
			puts "$size === $directory has previously been moved to individual --> removed folder from MSA"
			set directory_already_done_switch 1
			file delete -force $direct/MSA_split/$group/$size/[lindex $argv 0]/$directory
		}

		if {$directory_already_done_switch == 1} {continue} else {puts "$size === $directory is not yet DONE --> running analysis ..."}

		################################################################
		## Cd into directory and pick up current time for calcs below ##
		################################################################
		cd $direct/MSA_split/$group/$size/[lindex $argv 0]/$directory
		set current_time [clock seconds]

		##########################################################################################################
		## For each of the sizes, set the target number of independent ML reconstructions (the small should		##
		## typically use the -f a option with 100 BS replicates which then entails 20 ML trees - 1 per 5 BS).	##
		## The ML_tree_build_no can be modulated below 															##
		##########################################################################################################
		if {$size == "small" || $size == "large" || $size == "giant"} {
			set ML_tree_build_no 20
		} elseif {$size == "colossal"} {
			set ML_tree_build_no 10
		}
		
		#######################################################################
		## Check whether the BS file has been completed (small/large sizes). ##
		## If it has, rename it so we don't need to rerun the BS 			 ##
		#######################################################################
		set BS_done_switch 0
		set bootstrap_file [glob -nocomplain *bootstrap*]
		if {[llength $bootstrap_file] == 1} {
			set bootstrap_file [lindex $bootstrap_file 0]
			openfile $bootstrap_file
			set completed_BS [regexp -all {;} $data]
			if {$completed_BS == 100} {
				set BS_done_switch 1
				set temp_BS_name [string range $bootstrap_file 6 end]
				file rename $bootstrap_file $temp_BS_name
			}
		}

		##################################################################################################
		## Below we want to identify if there are any remnants from a previously aborted run - i.e.		##
		## some independently reconstructed ML trees, but the analysis did not have time to finish.		##
		## If this is the case - this could be for two reasons:											##
		##			(1) There was simply too little time, but given the full walltime it could finish	##
		## 			(2) The full walltime would not sufficient to complete the requested number of 		##
		## 			   of independent ML trees 															##
		## So, we need to see how long each ML trees, and see if we can fit the reconstruction into 	##
		## the walltime. If not, we try to incrementally decrease the requested number of ML trees 		##
		## to fit into the walltime with a minimum of AT LEAST 10 ML trees (same as for colossal). 		##
		##################################################################################################

		set time_per_tree 0

		## Find any incompleted Raxml intermediate ML tree files ##
		set ML_result_trees [glob -nocomplain RAxML_result*]
		
		## As long as there are at least two files, we can figure out how long it takes to build one ML tree. We cannot use the last result file as it is continuously updated and its timestamp may not represent the completed tree. We use the penultimate completed results tree - hence we need at least two ##

		## If possible, recalculate the time per tree - if the parsimony RUN.0 file exists##
		
		set parsimony_0_file [glob -nocomplain *parsimony*RUN.0]

		if {[llength $ML_result_trees] > 2 && [llength $parsimony_0_file] != 0} {
			set number_ML_completed [expr [llength [split [exec ls -lt | grep result] \n]] -1]
			set first_parsimony_file [lindex [split [exec ls -lt | grep parsimony | grep RUN.0]] end]; set first_parsimony_file_t [exec stat -c%Y $first_parsimony_file]
			set last_tree_file [lindex [split [lindex [split [exec ls -lt | grep result] \n] 1] " "] end]; set last_tree_file_t [exec stat -c%Y $last_tree_file]
			set total_time_taken [expr ([exec stat -c%Y $last_tree_file] - [exec stat -c%Y $first_parsimony_file]) / double(3600)]
			set time_per_tree [expr roundto([expr $total_time_taken / $number_ML_completed], 2)]
			# Write an output log file to record the average time-per-ML tree #
			set out [open ML_tree_time_log.txt w]
			puts $out "Average time per ML_tree for $directory based on $number_ML_completed completed ML trees:\t$time_per_tree"
			close $out
		## If the calculation has already been done before, and there are no results files to recalculate --> read the old calculated value ##
		} elseif {[file exists ML_tree_time_log.txt] == 1} {
			openfile ML_tree_time_log.txt
			set data [string trim $data]
			if {[string length $data] != 0} {
				set time_per_tree [lindex [split $data \t] end]
			} else {
				file delete ML_tree_time_log.txt
			}
		}
		## If we have a calculated time per ML tree we can adjust the number of trees to be calculated ##
		if {$time_per_tree != 0} {
			## For large or giant trees (or rerun of a small tree which already had a bootstrap file) we wanted 20 distinct ML trees that could be used to find the best one ##
			if {$size == "small" || $size == "large" || $size == "giant"} {
				## If building 20 trees would exceed the soft walltime limit (- 3 hours for variation in run time and final optimisation), we want to reduce the number of trees to be built ##
				if {[expr $time_per_tree * 20] > [expr $walltime_limit - 3]} {
					## If we can build 15 trees in the super-soft walltime, set the number to 15 ##
					if {[expr $time_per_tree * 15] < [expr $walltime_limit - 3]} {
						set ML_tree_build_no 15
						puts "This alignment: $group/$size/$directory takes $time_per_tree hours to construct a single ML tree. Number ML trees to be built: $ML_tree_build_no"
					## If we can build 10 trees in the super-soft walltime, set the number to 10 ##
					} elseif {[expr $time_per_tree * 10] < [expr $walltime_limit - 3]} {
						set ML_tree_build_no 10
						puts "This alignment: $group/$size/$directory takes $time_per_tree hours to construct a single ML tree. Number ML trees to be built: $ML_tree_build_no"
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
			}
		}

		##########################################################################
		## Once the ML tree time has been calculated (above), or if this was	##
		## an abortive -f a run or otherwise, delete the old raxml-specific 	##
		## files and rename the BS file (if exists) back						##
		##########################################################################
		set tree_build_files [glob -nocomplain RAxML*]
		foreach file $tree_build_files {file delete $file}

		## Once the other raxml-related files have been deleted, rename the BS file back to the original nomenclature ##
		if {$BS_done_switch == 1} {file rename $temp_BS_name $bootstrap_file}

		######################################################################
		## Calculated projected run time & exit early if exceeding walltime ##
		######################################################################

		## If we know time per tree and no. trees to be built - we can calculate projected run time ##
		if {$time_per_tree != 0} {
			# In hours
			set projected_run_time_hrs [expr $ML_tree_build_no * $time_per_tree]
			puts "The projected run-time for this tree (in hrs) based on previous runs:\t$projected_run_time_hrs"
			# In seconds
			set projected_run_time_secs [expr $projected_run_time_hrs * 3600]
			set projected_finish_time [expr $current_time + $projected_run_time_secs]
		## If we don't - use the shit method by simply multiplying the previous run time by a factor of 2 ##
		} else {
			set previous_run_time [expr $current_run_finish - $current_run_start]
			set projected_run_time_secs [expr round($previous_run_time * 2)]
			set projected_finish_time [expr $current_time + $projected_run_time_secs]
		}

		# If there is insufficient time left for the next time, exit #
		if {$projected_finish_time > $soft_end_time} {
			puts "DONE: $projected_finish_time > $soft_end_time ==== $current_time"
			exit 2
		}
		set current_run_start [clock seconds]

		#########################################
		## Process the alignment and key files ##
		#########################################
		set align [lindex [glob *.fas] 0]
		regsub {.fas} $align {} out_number
		set key [lindex [glob KEY*] 0]
		openfile $key; set key $data

		##########################################################################
		## Run the raxml tree reconstructions differentially dependent on $size ##
		##########################################################################
		# For small, we want 100 BS replicates and an ML tree constructed from at least 20 ML attempts #
		if {$size == "small"} {
			# If the BS replicates have not been completed, perform full analysis #
			if {$BS_done_switch == 0} {
				catch {exec raxml -f a -s $align -n $out_number\.txt -m PROTCATAUTO -p 1234 -x 10001 -N 100 -T $cores}
			# If they have, perform just the ML search (WITH CONVERGENCE BECAUSE SOMETIMES THEY ARE JUST TOO BIG/GAPPY) and then combine the BS and best ML trees #
			} else {
				catch {exec raxml -D -s $align -n $out_number\.txt -m PROTCATAUTO -p 1234 -N $ML_tree_build_no -T $cores}
				set ML_best [lindex [glob *bestTree*] 0]
				set BS [lindex [glob *bootstrap*] 0]
				# Assemble the final tree from the ML best tree and the BS replicates #
				catch {exec raxml -f b -t $ML_best -z $BS -m PROTCATAUTO -n STITCH_$out_number\.txt -T $cores}
			}
		# For large we want the same as for "small", but to use the convergence criterion #
		} elseif {$size == "large"} {
			if {$BS_done_switch == 0} {
				# Perform 100 rapid BS #
				catch {exec raxml -x 10001 -p 1234 -s $align -n BS_$out_number\.txt -m PROTCATAUTO -N 100 -T $cores}
			}
			# Perform a 20 ML search with convergence criterion = ON #
			catch {exec raxml -D -s $align -n ML_$out_number\.txt -m PROTCATAUTO -p 1234 -N $ML_tree_build_no -T $cores}
			set ML_best [lindex [glob *bestTree*] 0]
			set BS [lindex [glob *bootstrap*] 0]
			# Assemble the final tree from the ML best tree and the BS replicates #
			catch {exec raxml -f b -t $ML_best -z $BS -m PROTCATAUTO -n STITCH_$out_number\.txt -T $cores}
		# For the giant and colossal categories, we only do ML best trees based on 20 and 10 ML trees respectively #
		} else {
			if {$size == "giant"} {
				catch {exec raxml -D -s $align -n $out_number\.txt -m PROTCATAUTO -p 1234 -N $ML_tree_build_no -T $cores}
			} elseif {$size == "colossal"} {
				catch {exec raxml -D -s $align -n $out_number\.txt -m PROTCATAUTO -p 1234 -N $ML_tree_build_no -T $cores}
			}
		}

		##############################################################################################
		## If a file with BS values is available (small/large) - use that to rename the tip labels.	##
		## Otherwise use the best ML tree. 															##
		##############################################################################################
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
			file copy -force $out_number\_tree.txt $direct/Trees/$size/

			cd $direct/MSA_split/$group/$size/[lindex $argv 0]
			file rename -force $directory $direct/Tree_build/$size/
		# If there is no "best tree" escape the script #
		} else {
			puts "There were no completed tree files ... Walltime exceeded?"
			exit 2 
		}
		##################################################################
		## Set the finish time for the current directory reconstruction ##
		##################################################################
		set current_run_finish [clock seconds]
	}
}