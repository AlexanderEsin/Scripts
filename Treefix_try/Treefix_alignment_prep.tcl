## Prepare the alignment files for analysis with treefix-DTL ##

source ~/Dropbox/Scripts/General_utils.tcl

###########################################################################

set direct /users/aesin/desktop/Test_dtl_methods/Ranger/Treefix_test

set trees_to_test_dir $direct/Final_BS
set pre_dir $direct/Alignment_pre
set post_dir $direct/Alignment_post
set new_label_key_dir $direct/Tip_keys

cd $trees_to_test_dir
set trees [glob *_tree.txt]

foreach tree $trees {
	set number [string range $tree 0 [string first \_ $tree]-1]

	## Reformat the alignment from phylip-interleaved into fasta format ##
	cd $pre_dir
	set alignment "$number\.fas"
	catch {exec seqret $alignment $number\_fasta.fasta}

	## Open the two keys. 1: from the k-val format to the binomial_class_proteinID. 2: from binomial_class_proteinID to binomial_number_taxa
	set k_val_key_data [split [string trim [openfile $pre_dir/KEY_$number\.txt]] \n]
	set tip_key_data [split [string trim [openfile $new_label_key_dir/TIP_KEY_$number\.tsv]] \n]

	## Open the fasta-formatted alignment and split into individual sequences ##
	set alignment_seqs [split_genes [string trim [openfile $pre_dir/$number\_fasta.fasta]]]
	puts "Number of sequences: [llength $alignment_seqs]"

	set new_align {}

	foreach seq $alignment_seqs {
		set seq_name [string trim [string range $seq 1 [string first \n $seq]]]
		
		set long_tip_values [split [string trim [lsearch -inline -glob $k_val_key_data $seq_name\t*]] \t]
		set long_tip_name "[lindex $long_tip_values 1]\{[lindex $long_tip_values 2]\}\{[lindex $long_tip_values 3]\}"

		regsub -all {__} $long_tip_name "_" long_tip_name

		set final_tip_values [split [string trim [lsearch -inline -glob $tip_key_data $long_tip_name*]] \t]
		set final_tip_name [lindex $final_tip_values 1]

		set new_seq ">$final_tip_name\n[string trim [string range $seq [string first \n $seq] end]]"
		lappend new_align $new_seq
	}

	set out [open $post_dir/Align_$number\.fasta w]
	puts $out [join $new_align \n]
	close $out
}