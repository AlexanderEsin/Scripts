#!/usr/bin/tclsh

set direct /home/ade110/Work/Bacillus/Prefilter/Bacillus_withOut

set tmpdir $::env(TMPDIR)

set group_dir	$direct/BLAST

set all_in_dir	$group_dir/IN_all_split
set all_db_dir	$group_dir/DB_all

## Get subfolder
set fold_no 	[lindex $argv 0]
set fold_name	"dir_$fold_no"

## Make output
set out_main	$group_dir/OUT_all_split2
set out_sub		$out_main/$fold_name
file mkdir 		$out_sub


######################

set proteome_l	[glob $all_in_dir/$fold_name/*.fasta]
set all_db_l	[glob $all_db_dir/*.pin]

foreach proteome $proteome_l {

	## Reduce to just acc_ass
	set f_query_file	[file tail $proteome]
	set f_query_base	[string range $f_query_file 0 end-15]

	## Each proteome is blasted against all other proteomes (dbs) in the same group
	foreach db $all_db_l {

		## Reduce to just acc_ass
		set f_db_file	[file tail $db]
		set f_db_base 	[string range $f_db_file 0 end-4]

		if {$f_query_base != $f_db_base} {
			exec -ignorestderr blastp -query $proteome -db "$all_db_dir/$f_db_base" -out "$out_sub/$f_query_base\&$f_db_base\.tsv" -evalue 1e-10 -outfmt 6 -max_target_seqs 1 -max_hsps 1 -seg "yes" -soft_masking true -num_threads 1
		}
	}
}