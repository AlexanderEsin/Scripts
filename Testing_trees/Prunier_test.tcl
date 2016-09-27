source ~/Dropbox/Scripts/General_utils.tcl

set direct ~/desktop/Prunier_test
# file mkdir $direct/Relabelled
# file mkdir $direct/Pruned

cd $direct
# set trees [glob *_tree.txt]

# foreach tree $trees {
# 	set tree_num [string range $tree 0 [string first \_ $tree]-1]

# 	set tree_data [string trim [openfile $tree]]

# 	set pata {(?:\s+\{)+?.+?\}+?}

# 	regsub -all $pata $tree_data {} new_tree_data

# 	set out [open $direct/Relabelled_$tree_num\.txt w]
# 	puts $out [string trim $new_tree_data]
# 	close $out
# }

set input_tree_file [lindex [glob Relabelled_*] 0]
set input_tree_data [string trim [openfile $input_tree_file]]

set reference_tree_file [lindex [glob Astral* ] 0]
set reference_tree_data [string trim [openfile $reference_tree_file]]

set alns [glob *fas]
set new_aln_file {}
set binomial_l {}

foreach aln $alns {
	set number [regexp -inline {[0-9]+} $aln]
	catch {exec seqret $aln -out $number\_fasta.fas}
	after 1000

	set proteins [split_genes [string trim [openfile $number\_fasta.fas]]]
	set key_data [split [string trim [openfile KEY_$number\.txt]] \n]

	foreach protein $proteins {
		set header [string trim [string range $protein 0 [string first \n $protein]]]
		set k_val [string range $header 1 end]
		set key_hit [lsearch -inline -glob $key_data $k_val*]
		set key_entries [split $key_hit \t]
		set binomial [lindex $key_entries 1]

		set new_entry ">$binomial\n[string range [string trim $protein] [string first \n $protein]+1 end]"
		lappend new_aln_file $new_entry
	}

	set new_aln_file_str [join $new_aln_file \n]
	set out [open $aln\.txt w]; puts $out $new_aln_file_str; close $out

	foreach entry $new_aln_file {

		set header_binomial [string range [string trim $entry] 1 [string first \n $entry]-1]

		lappend binomial_l $header_binomial

		if {[regexp $header_binomial $reference_tree_data] == 0} {
			puts "Not found in reference: $header_binomial"
		}
		if {[regexp $header_binomial $input_tree_data] == 0} {
			puts "Not found in gene tree: $header_binomial"
		}
	}

	set duplicates [dups $binomial_l]

	file delete $number\_fasta.fas
}

set aln_file $aln\.txt


chan configure stdout -buffering none
exec prunier input.tree.file=$reference_tree_file aln.file=$aln_file sequence.type=protein raxml.path=/usr/bin/raxml genetree.file=$input_tree_file raxml.nb_proc=4 aln.type=FASTA