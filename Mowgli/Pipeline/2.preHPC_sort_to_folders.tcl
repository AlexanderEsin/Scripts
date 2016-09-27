#!/usr/local/bin/tclsh
## Run mowgli on the HPC ##

## For desktop ##
source ~/Dropbox/Scripts/General_utils.tcl
set direct /users/aesin/Desktop/Mowgli

# ## For HPC ##
# source /home/ade110/Scripts/General_utils.tcl
# set direct /csc/rawdata/Warnecke/ADE/Mowgli

set t_values {3 4 5 6}

###########################################################################

set input_directory Tips_relabelled_BS

foreach t_value $t_values {
	set output_directory HPC_mowgli_input_BS_$t_value
	
	set input_path $direct/$input_directory
	set output_path $direct/$output_directory

	file mkdir $output_path

	## Get list of all tree files in the input directory ##
	cd $input_path
	set trees_to_test [lsort -dictionary [glob *txt]]

	set folder_counter 1

	## Foreach tree, make a new directory to house the tree file ##
	foreach tree $trees_to_test {
		file mkdir $output_path/$folder_counter
		file copy $tree $output_path/$folder_counter

		incr folder_counter
		puts $folder_counter
	}

}


###########################################################################


