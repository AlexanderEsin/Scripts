## Run mowgli on the HPC ##

source ~/Scripts/General_utils.tcl

proc reconcile_success {} {
	upvar 1 tree_number tree_number
	upvar 1 temp_output_path temp_output_path
	upvar 1 success_path success_path

	file rename -force $temp_output_path/$tree_number $success_path/
	puts stdout "Reconciliation was successful: $tree_number"
	return
}

###########################################################################
## Set penalty cose ##
set transfer_cost 8

###########################################################################
## Set species tree location and names of directories ##
set direct /csc/rawdata/Warnecke/ADE/Mowgli
set species_tree_file $direct/Species_tree/Ultrametric_species_tree_tipfix.txt

set input_directory HPC_mowgli_input_BS_$transfer_cost
set failed_directory HPC_input_BS_FAILED_$transfer_cost
set success_directory Mow_test_out_t$transfer_cost

###########################################################################
## Set paths and make necessary directories ##
set input_path $direct/$input_directory
set failed_path $direct/$failed_directory
set temp_output_path $direct/Mowgli_outputs/Temp
set success_path $direct/Mowgli_outputs/$success_directory

file mkdir $failed_path
file mkdir $temp_output_path
file mkdir $success_path

####################################
## Navigate to input and get the tree ##
set folder_number [lindex $argv 0]
cd $input_path/$folder_number
set input_tree [lindex [glob *txt] 0]

set tree_number [regexp -inline {[0-9]+} $input_tree]

####################################
## Prepare output and run Mowgli ##
puts stdout "Testing tree: $tree_number"
file mkdir $temp_output_path/$tree_number

# Model 1 (-M 1 is default)
# m (number of MPRs computed - default = "1")
# n (gene tree modification by NNI. 1 = "ON")
# T (threshold for BS support for NNI modification)
set success [catch {exec mowgli -s $species_tree_file -g $input_tree -v -n 1 -T 80 -o $temp_output_path/$tree_number/ -t $transfer_cost -d 2 -l 1} errvar]

if {[string trim $errvar] eq "Exception: NNI() ----> invalid NNI"} {
	## Linux Mowgli can't deal with trying to do NNI when NNI cannot be done - so rerun without NNI algorithm ##
	puts stdout "Re-running mowgli with NNI off..."
	set success [catch {exec mowgli -s $species_tree_file -g $input_tree -v -n 0 -T 80 -o $temp_output_path/$tree_number/ -t $transfer_cost -d 2 -l 1} errvar]
}

if {$success == 0} {
	reconcile_success
} elseif {[string trim $errvar] eq "NO BETTER GENE TREE ARRANGEMENT FOUND."} {
	reconcile_success
} else {
	file delete -force $temp_output_path/$tree_number
	file copy -force $input_path/$folder_number/$input_tree $failed_path/
	puts stdout "Reconciliation failed: $tree_number\. Error code: $success\. Error message: $errvar"
}
	
####################################
