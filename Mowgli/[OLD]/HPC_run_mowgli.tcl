## Run mowgli on the HPC ##

source ~/Scripts/General_utils.tcl

set transfer_cost 3

set direct /csc/rawdata/Warnecke/ADE/Mowgli
set species_tree_file $direct/Species_tree/Ultrametric_species_tree_tipfix.txt

set input_directory HPC_mowgli_input_BS
set failed_directory HPC_input_BS_FAILED
set success_directory Mow_test_out_t$transfer_cost

###########################################################################

set input_path $direct/$input_directory
set failed_path $direct/$failed_directory
set temp_output_path $direct/Mowgli_outputs/Temp
set success_path $direct/Mowgli_outputs/Mow_test_out_t$transfer_cost

file mkdir $failed_path
file mkdir $temp_output_path
file mkdir $success_output_directory

####################################

set folder_number [lindex $argv 0]
cd $input_path/$folder_number
set input_tree [lindex [glob *txt] 0]

####################################

set tree_number [regexp -inline {[0-9]+} $input_tree]
puts stdout "Testing tree: $tree_number"
file mkdir $temp_output_path/$tree_number

## Model 1 (-M 1 is default)
## m (number of MPRs computed - default = "1")
## n (gene tree modification by NNI. 1 = "ON")
## T (threshold for BS support for NNI modification)
set success [catch {exec mowgli -s $species_tree_file -g $input_tree -v -n 1 -T 80 -o $temp_output_path/$tree_number/ -t $transfer_cost -d 2 -l 1}]

if {$success == 0} {
	file rename -force $temp_output_path/$tree_number $success_path/
	puts stdout "Reconciliation was successful"
} else {
	file delete -force $temp_output_path/$tree_number
	file copy -force $input_path/$folder_number/$input_tree $failed_path/
	puts stdout "Reconciliation failed. Error code: $success"
}
	
####################################