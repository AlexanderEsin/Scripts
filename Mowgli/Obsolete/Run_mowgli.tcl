#!/usr/local/bin/tclsh
## Trial run of mowgli ##


set parallel_threads 7
set process_number [lindex $argv 0]
set skip_factor $process_number

set transfer_cost 4


set direct /users/aesin/desktop/Test_dtl_methods/Mowgli
set species_tree_file $direct/Ultrametric_species_tree_tipfix.txt

set testing_directory $direct/Random_trees_21

###########################################################################

set input_tree_directory $testing_directory/Tips_relabelled
set output_directory $testing_directory/Mowgli_output_t$transfer_cost

file mkdir $output_directory

cd $input_tree_directory
set all_input_trees [glob *txt]

## Assignment of files for parellelised version ##
set input_trees {}
set files_parsed 0
while {$files_parsed < [llength $all_input_trees]} {
	lappend input_trees [lindex $all_input_trees $skip_factor]
	set skip_factor [expr $skip_factor + $parallel_threads]
	set files_parsed $skip_factor
}

puts "Process $process_number\n-->Trees to be tested by Mowgli: [llength $input_trees]"
#puts $input_trees

foreach tree $input_trees {
	set tree_number [regexp -inline {[0-9]+} $tree]
	puts "Process $process_number\n-->Testing tree: $tree_number"
	file mkdir $output_directory/$tree_number
	
	## Model 1 (-M 1 is default)
	## m (number of MPRs computed - default = "1")
	## n (gene tree modification by NNI. 1 = "ON")
	## T (threshold for BS support for NNI modification)
	catch {exec mowgli -s $species_tree_file -g $tree -v -n 1 -T 80 -o $output_directory/$tree_number/ -t $transfer_cost -d 2 -l 1}
}