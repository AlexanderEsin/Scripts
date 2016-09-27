##################
## Source procs ##
##################

source ~/Scripts/General_utils.tcl

###########################################################################

## The ival at which the clusteres were made ##
set ival 2.0
#set direct /users/aesin/desktop/Geo_v_all/$ival
set direct /scratch/ade110/Geo_v_all/$ival/Trees/Final_NBS

set cores 3
set incomplete_file_counters {}

###########################################################################
## Identify the folder in which the tree is and make a "Bootstraps" directory there ##
cd $direct
set tree_folder [lindex [lsort -dictionary [glob -type d * ]] [expr [lindex $argv 0] - 1]]
puts "Working on: $tree_folder"

file mkdir $direct/$tree_folder/Bootstraps

#####################################
## Check how many bootstraps we have already been made ##

cd $direct/$tree_folder/Bootstraps

set BS_output_files [lsort -dictionary [glob -nocomplain BS_output*]]
set incomplete_file_counters {}
set bs_files_counter 11

if {[llength $BS_output_files] == 0} {
	set bs_files_counter 1

} elseif {[llength $BS_output_files] < 10} {

	## Let's find all the files that do not exist ##
	for {set i 1} {$i < 11} {incr i} {
		if {[file exists BS_output_$i\.txt] == 1} {
			continue
		} else {
			lappend incomplete_file_counters $i
		}
	}
} else {
	## Check that each completed BS file has 10 bootstraps as expected ##
	foreach completed_file $BS_output_files {
		set bs_counter [regexp -all ";" [openfile $completed_file]]
		if {$bs_counter == 10} {
			puts "There are 10 bootstraps in $completed_file"
		} else {
			puts "There are fewer than 10 bootstraps $completed_file. Bootstrap number: $bs_counter"
			set file_number [string range $completed_file [string last \_ $completed_file]+1 [string first \. $completed_file]-1]
			lappend incomplete_file_counters $file_number
		}
	}
	puts $incomplete_file_counters
}

## Delete intermediate RAxML files ##
set raxml_files [glob -nocomplain RAxML*]
if {[llength $raxml_files] != 0} {foreach raxml_file $raxml_files {file delete $raxml_file}}

###########################################################################

while {$bs_files_counter <= 10} {
	cd $direct/$tree_folder
	set align [lindex [glob *fas] 0]

	set random_seed_x [RandomInt 1000 25000]
	set random_seed_p [RandomInt 1000 25000]

	catch {exec raxml -x $random_seed_x -p $random_seed_p -s $align -w $direct/$tree_folder/Bootstraps -n BS_$bs_files_counter\.txt -m PROTCATAUTO -N 10 -T $cores}

	cd $direct/$tree_folder/Bootstraps
	file rename RAxML_bootstrap.BS_$bs_files_counter\.txt BS_output_$bs_files_counter\.txt

	set raxml_files [glob RAxML*]
	foreach raxml_file $raxml_files {file delete $raxml_file}

	incr bs_files_counter
}

###########################################################################

if {[llength $incomplete_file_counters] > 0} {
	foreach counter $incomplete_file_counters {

		## Delete the file with the incomplete bootstraps ##
		file delete BS_output_$counter\.txt

		## Find the alignment ##
		cd $direct/$tree_folder
		set align [lindex [glob *fas] 0]

		## Get random seeds ##
		set random_seed_x [RandomInt 1000 25000]
		set random_seed_p [RandomInt 1000 25000]

		## Run bootstraps ##
		catch {exec raxml -x $random_seed_x -p $random_seed_p -s $align -w $direct/$tree_folder/Bootstraps -n BS_$counter\.txt -m PROTCATAUTO -N 10 -T $cores}

		cd $direct/$tree_folder/Bootstraps
		file rename RAxML_bootstrap.BS_$counter\.txt BS_output_$counter\.txt

		set raxml_files [glob RAxML*]
		foreach raxml_file $raxml_files {file delete $raxml_file}
	}
}