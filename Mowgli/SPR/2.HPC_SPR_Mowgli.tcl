## Run mowgli on the HPC ##

source ~/Scripts/General_utils.tcl

proc reconcile_success {} {
	upvar 1 tree_number tree_number
	upvar 1 temp_output_path temp_output_path
	upvar 1 success_path success_path
	upvar 1 spr_num spr_num

	file rename -force $temp_output_path/$tree_number/$spr_num $success_path/
	puts stdout "Reconciliation was successful: $tree_number SPR $spr_num"
	return
}

###########################################################################
## Set penalty cose ##
set transfer_cost 5

###########################################################################
## Set species tree location and names of directories ##
set direct /csc/rawdata/Warnecke/ADE/Mowgli
set species_tree_file $direct/Species_tree/Ultrametric_species_tree_tipfix.txt

set input_directory SPR_mowgli_input
set failed_directory SPR_mow_FAILED_$transfer_cost
set success_directory SPR_mow_out_t$transfer_cost

###########################################################################
## Set paths and make necessary directories ##
set input_path $direct/$input_directory
set failed_path $direct/$failed_directory
set temp_output_path $direct/Mowgli_outputs/Temp
set success_path $direct/Mowgli_outputs/$success_directory

file mkdir $failed_path
file mkdir $temp_output_path

####################################
## Navigate to input and get the tree ##
set iteration_number [expr [lindex $argv 0] - 1]
#set iteration_number 1
cd $input_path
set tree_number [lindex [glob -type d *] $iteration_number]

set failed_path $failed_path/$tree_number
set success_path $success_path/$tree_number
file mkdir $success_path

cd $input_path/$tree_number
set SPR_trees [lsort -dictionary [glob *txt]]

foreach spr_tree $SPR_trees {
	set spr_num [regexp -inline {[0-9]+} $spr_tree]

	puts stdout "Testing tree: $tree_number SPR: $spr_num"

	file mkdir $temp_output_path/$tree_number/$spr_num

	if {[file exists $temp_output_path/$tree_number/$spr_num/costs.mpr] == 1} {
		if {[regexp "Costs" [openfile $temp_output_path/$tree_number/$spr_num/costs.mpr]] > 0} {
			reconcile_success
			continue
		}
	}

	set success [catch {exec mowgli -s $species_tree_file -g $spr_tree -v -n 1 -T 80 -o $temp_output_path/$tree_number/$spr_num -t $transfer_cost -d 2 -l 1} errvar]

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
		file mkdir $failed_path

		file delete -force $temp_output_path/$tree_number/$spr_num
		file copy -force $input_path/$tree_number/$spr_tree $failed_path/
		puts stdout "Reconciliation failed: $tree_number SPR: $spr_num\. Error code: $success\. Error message: $errvar"
	}

}

