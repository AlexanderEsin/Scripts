#!/usr/local/bin/tclsh
#############
### Procs ###
#############

source ~/Dropbox/Scripts/General_utils.tcl

###########################################################################

set direct /users/aesin/Desktop/Geo_analysis/Geo_v_all/2.0/Trees/NBS_to_BS_HPC

###########################################################################

cd $direct
set new_BS_groups [glob -type d *]

puts "Folders to be processed [llength $new_BS_groups]"
set master_BS_made 0

## Regexp patterns ##
set pata {k+?[0-9]+?:+?}
set patb {\t+?.+?\n+?}

#####################################
#set new_BS_groups [lrange $new_BS_groups 0 0]
foreach new_BS_group $new_BS_groups {
	cd $direct/$new_BS_group/Bootstraps

	#####################################
	## Each file contains 10 BS replicates ##
	set BS10_files [glob -nocomplain -type f *]
	set number_of_files [llength $BS10_files]
	if {$number_of_files != 10} {
		puts "For group $new_BS_group the number of BS repeats is less than 100. There are $number_of_files files done"
		continue
	}

	#####################################
	## Make a master bootstrap file with 100 BS replicates ##
	set master_BS_file {}
	foreach BS10_file $BS10_files {
		set 10_bootstraps [string trim [openfile $BS10_file]]
		lappend master_BS_file $10_bootstraps
	}

	set master_BS_file [join $master_BS_file \n]
	set out [open $direct/$new_BS_group/master_BS_file.txt w]
	puts $out $master_BS_file
	close $out


	#####################################
	## Stitch the bootstraps onto the best_tree ##
	cd $direct/$new_BS_group

	set ML_best [lindex [glob *bestTree*] 0]
	catch {exec raxml -f b -t $ML_best -z master_BS_file.txt -m PROTCATAUTO -n STITCH_$new_BS_group\.txt -T 4}

	#####################################
	## Rename the old translated tip tree file to mark is at NBS ##
	set old_NBS_file "$new_BS_group\_tree.txt"
	if {[file exists $old_NBS_file] == 1} {
		file rename $old_NBS_file "NBS_$old_NBS_file"
	}

	#####################################
	## Translate the tip labels of the stitched output tree ##
	set tree_data [openfile RAxML_bipartitions.STITCH_$new_BS_group\.txt]
	set branchids [regexp -all -inline $pata $tree_data]

	set key [lindex [glob KEY*] 0]
	set key [openfile $key]

	foreach k_val $branchids {
		regsub {:} $k_val {} k_val_trunc
		regexp -line $k_val_trunc$patb $key hit
		set binomial [lindex [split [string trim $hit] \t] 1]
		set class [lindex [split [string trim $hit] \t] 2]
		set id [lindex [split [string trim $hit] \t] 3]
		# If you use k_val_trunc below, you may substitute k17 when you want to find k1 - so searching for k1: is more accurate #
		regsub $k_val $tree_data "$binomial \{$class\} \{$id\}:" tree_data
	}

	set out [open BS_$new_BS_group\_tree.txt w]; puts $out $tree_data; close $out

	incr master_BS_made
}

puts "Number of master BS files created = $master_BS_made / [llength $new_BS_groups]"