#!/usr/bin/tclsh

## Run mowgli on the HPC ##

source /home/ade110/Scripts/General_utils.tcl

proc reconcile_success {} {
	upvar 1 tree_number tree_number
	upvar 1 temp_output_path temp_output_path
	upvar 1 success_path success_path

	file rename -force $temp_output_path/$tree_number $success_path/
	puts stdout "Reconciliation was successful: $tree_number"
	return
}

## Set penalty cost
set transfer_cost 	[lindex $argv 1]

## Set species tree location and names of directories
set direct			/scratch/ade110/FastTree
set species_tree	$direct/Species_tree/Ultrametric_species_tree_tipfix.txt

set input_dir		$direct/Mowgli_input/Input_$transfer_cost
set failed_dir		$direct/Mowgli_failed/Failed_$transfer_cost
set success_dir		$direct/Mowgli_output/Output_$transfer_cost
set temp_dir		$direct/Mowgli_output/Temp

## Make necessary directories
file mkdir $failed_dir $success_dir $temp_dir


####################################
## Navigate to input and get the tree ##
set folder_number	[lindex $argv 0]

set tree_file		[lindex [glob $input_dir/$folder_number/*txt] 0]
set file_name		[file tail $tree_file]
set tree_number		[string range $file_name [string last \_ $file_name]+1 end-4]


####################################
## Prepare output and run Mowgli ##
puts stdout			"Testing tree: $tree_number"
set output_dir		$success_dir/$tree_number
file mkdir 			$output_dir

# Model 1 (-M 1 is default)
# m (number of MPRs computed - default = "1")
# n (gene tree modification by NNI. 1 = "ON")
# T (threshold for BS support for NNI modification)
set success [catch {exec mowgli -s $species_tree -g $tree_file -v -n 1 -T 80 -o $output_dir/ -t $transfer_cost -d 2 -l 1} errvar]
puts stdout $success
puts stdout $errvar

# if {[string trim $errvar] eq "Exception: NNI() ----> invalid NNI"} {
# 	## Linux Mowgli can't deal with trying to do NNI when NNI cannot be done - so rerun without NNI algorithm ##
# 	puts stdout "Re-running mowgli with NNI off..."
# 	puts stdout "Re-running mowgli with NNI off..."
# 	set success [catch {exec mowgli -s $species_tree_file -g $input_tree -v -n 0 -T 80 -o $temp_output_path/$tree_number/ -t $transfer_cost -d 2 -l 1} errvar]
# }

# if {$success == 0} {
# 	reconcile_success
# } elseif {[string trim $errvar] eq "NO BETTER GENE TREE ARRANGEMENT FOUND."} {
# 	reconcile_success
# } else {
# 	file delete -force $temp_output_path/$tree_number
# 	file copy -force $input_path/$folder_number/$input_tree $failed_path/
# 	puts stdout "Reconciliation failed: $tree_number\. Error code: $success\. Error message: $errvar"
# }
	
####################################
