## Run mowgli on the HPC ##
source /home/ade110/Scripts/General_utils.tcl

## Set penalty cost
set transfer_cost 	[lindex $argv 0]
#set transfer_cost 	3


## Set species tree location and names of directories
set direct			/home/ade110/Work/Bacillus
set mowgli_dir		$direct/Mowgli
set species_tree	$direct/Astral_speciesTree/Astral_speciesTree_ultraChecked.txt

set input_dir		$mowgli_dir/GeneTree_input/
set input_list		[lreverse [lsort -dictionary [glob -type d $input_dir/*]]]
set success_dir		$mowgli_dir/Mowgli_output/Output_$transfer_cost

## Make output directories
file mkdir $success_dir

####################################
## Navigate to input and get the tree ##
set folder_index	[expr [lindex $argv 1] - 1]
set input_folder	[lindex $input_list $folder_index]

set tree_file		[lindex [glob $input_folder/*FT_relabelled.txt] 0]
set file_name		[file tail $tree_file]
set tree_number		[string range $file_name 0 [string first \_ $file_name]-1]

####################################
## Prepare output and run Mowgli ##
puts stdout			"Testing tree: $tree_number"
set output_dir		$success_dir/$tree_number
file mkdir 			$output_dir

# If .ok file exists, reconciliation already succeeded
if {[file exists $output_dir/reconciled.ok] == 1} {
	puts stdout "\t>>> Already reconciled!"
	flush stdout
	exit
}

puts stdout "\t>>> Attempting to reconcile..."
set log	[open $output_dir/reconciled.log w]

# Model 1 (-M 1 is default)
# m (number of MPRs computed - default = "1")
# n (gene tree modification by NNI. 1 = "ON")
# T (threshold for BS support for NNI modification)
set success [catch {exec mowgli -s $species_tree -g $tree_file -v -n 1 -T 80 -o $output_dir -t $transfer_cost -d 2 -l 1} errvar]

puts $log $errvar

if {[file size $output_dir/Fullreconciliation.mpr] == 0} {
	puts $log "Reconciliation failed: $tree_number"
} else {
	puts $log "Reconciliation succeded: $tree_number"
	set ok [open $output_dir/reconciled.ok w]
	close $ok
}

flush stdout
close $log



# if {[file exists $output_dir/Fullreconciliation.mpr] == 1} {
# 	puts stdout "Reconciliation already attempted..."
# 	if {[file size $output_dir/Fullreconciliation.mpr] > 0} {
# 		puts stdout "\t... and was successful. Skipping"
# 		exit
# 	} else {
# 		puts stdout "\.... and was unsuccessful. Attempting again."
# 	}
# } else {
# 	puts stdout "Reconciliation not yet attempted..."
# }



# puts $errvar

# if {[file size $output_dir/Fullreconciliation.mpr] == 0} {
# 	puts stderr "Reconciliation failed: $tree_number"
# } else {
# 	puts stdout "Reconciliation success: $tree_number"
# }