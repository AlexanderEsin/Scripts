## Run mowgli on the HPC ##
source /mnt/storage/home/aesin/Scripts/General_utils.tcl

## Set penalty cost
set transfer_cost 	[lindex $argv 0]

## Set species tree location and names of directories
set direct			/mnt/storage/rawdata/Warnecke_Rawdata/ADE/Work/Staph
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
# puts stdout			"Testing tree: $tree_number"
set output_dir		$success_dir/$tree_number
file mkdir 			$output_dir

# If .ok file exists, reconciliation already succeeded
if {[file exists $output_dir/reconciled.ok] == 1} {
	puts stdout "Testing tree: $tree_number\n\t>>> Already reconciled!"
	flush stdout
	exit
}

puts stdout "Testing tree: $tree_number\n\t>>> Attempting to reconcile..."
flush stdout
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