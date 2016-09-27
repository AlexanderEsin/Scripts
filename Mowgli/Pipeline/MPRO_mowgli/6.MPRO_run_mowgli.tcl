## Run mowgli on the Mac Pro ##

source ~/Dropbox/Scripts/General_utils.tcl

set transfer_cost 6

set direct /users/aesin/desktop/Mowgli
set species_tree_file $direct/Species_tree/Ultrametric_species_tree_tipfix.txt

# set input_directory MPRO_failed/MPRO_input_BS_FAILED_$transfer_cost
set input_directory MPRO_input/MPRO_input_BS_$transfer_cost
set failed_directory MPRO_failed/BS_FAILED_$transfer_cost
set success_directory Mow_test_out_t$transfer_cost

###########################################################################

set input_path $direct/$input_directory
set failed_path $direct/$failed_directory
set temp_output_path $direct/Mowgli_outputs/Temp
set success_path $direct/Mowgli_outputs/$success_directory

file mkdir $failed_path
file mkdir $temp_output_path
file mkdir $success_path

####################################

set folder_number [lindex $argv 0]
cd $input_path/$folder_number

set input_trees [reverse_dict [lsort -dictionary [glob *txt]]]

####################################

foreach input_tree $input_trees {
	set tree_number [regexp -inline {[0-9]+} $input_tree]
	puts stdout "Testing tree: $tree_number"

	if {[file exists $success_path/$tree_number] == 1 || [file exists $failed_path/$input_tree] == 1} {
		puts stdout "Mowgli already completed (or failed to complete) the reconciliation, skipping ..."
		continue
	}

	file mkdir $temp_output_path/$tree_number

	## Model 1 (-M 1 is default)
	## m (number of MPRs computed - default = "1")
	## n (gene tree modification by NNI. 1 = "ON")
	## T (threshold for BS support for NNI modification)
	set success [catch {exec mowgli -s $species_tree_file -g $input_tree -v -n 1 -T 80 -o $temp_output_path/$tree_number/ -t $transfer_cost -d 2 -l 1} errvar]

	if {$success == 0} {
		file rename -force $temp_output_path/$tree_number $success_path/
		puts stdout "Reconciliation was successful: $tree_number"
	} else {
		## The 64-bit Mac binary throws this error when NNI can't be done, but it still completes the reconciliation ##
		if {[string trim $errvar] eq "NO BETTER GENE TREE ARRANGEMENT FOUND."} {
			file rename -force $temp_output_path/$tree_number $success_path/
			puts stdout "Reconciliation was successful: $tree_number"
		} else {
			file delete -force $temp_output_path/$tree_number
			file copy -force $input_path/$folder_number/$input_tree $failed_path/
			puts stdout "Reconciliation failed: $tree_number\. Error code: $success\. Error message: $errvar"
		}
	}	
}

####################################