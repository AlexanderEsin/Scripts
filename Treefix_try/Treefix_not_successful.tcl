source ~/Dropbox/Scripts/General_utils.tcl

###########################################################################

set direct /users/aesin/desktop/Test_dtl_methods/Ranger/Treefix_test

set trees_to_test_dir $direct/Tips_relabelled
set align_dir $direct/Alignment_post
set spec_map_dir $direct/Species_maps
set treefix_in_dir $direct/Treefix_input

###########################################################################

set align_extension ".fasta"

cd $trees_to_test_dir
set input_trees [glob *txt]

foreach tree $input_trees {
	set number [regexp -inline {[0-9]+} $tree]

	file copy -force $tree $treefix_in_dir/
	file rename -force $treefix_in_dir/$tree $treefix_in_dir/$number\.tree

	set align_file $align_dir/Align_$number$align_extension
	file copy -force $align_file $treefix_in_dir/
	file rename -force $treefix_in_dir/Align_$number$align_extension $treefix_in_dir/$number\.align


	set smap $spec_map_dir/SMAP_$number\.smap
	set species_tree $direct/Species_tree_gaps.txt

	cd $treefix_in_dir
	puts [exec pwd]
	chan configure stdout -buffering none
	puts "$number\.tree"
	

	# exec treefix_compute --type likelihood -m treefix.models.raxmlmodel.RAxMLModel -U .tree -e "-m PROTCATWAG -e 2.0" $number\.tree

	exec /bin/bash -c "treefix_compute --type cost -m treefix.models.rangerdtlmodel.DTLModel -s $species_tree -S $smap $number\.tree"

	# exec /bin/bash -c "treefix_compute --type cost -m treefix.models.rangerdtlmodel.DTLModel -s $species_tree -S $smap -o .tree -niter=100 $number\.tree"

	#exec treefixDTL -s $species_tree -S $smap -e "-m PROTCATWAG -e 2.0" -U .tree --niter=100 -E "--tmp ./Temp" -l log.txt -V3 $number\.tree
}

