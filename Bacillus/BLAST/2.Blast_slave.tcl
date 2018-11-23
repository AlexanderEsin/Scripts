#!/usr/bin/tclsh

set direct /home/ade110/Work/Bacillus/BLAST

set all_in_dir $direct/IN_all_split
set all_db_dir $direct/DB_all

set grp_in_dir $direct/IN_grp
set grp_db_dir $direct/DB_grp

## Get subfolder
set fold_no 	[lindex $argv 0]
set fold_name	"dir_$fold_no"


## Make output
set out_main	$direct/OUT_all_split
set out_sub		$out_main/$fold_name
file mkdir 		$out_sub


######################

set proteome_l	[glob $all_in_dir/$fold_name/*.fasta]
set grp_db_l	[glob $grp_db_dir/*.pin]
set blast_run	false


foreach proteome $proteome_l {

	## Reduce to just acc_ass
	set f_query_file	[file tail $proteome]
	set f_query_base	[string range $f_query_file 0 end-15]

	## Each proteome is blasted against the 25 group databases
	foreach db $grp_db_l {

		## Reduce to just acc_ass
		set f_db_file	[file tail $db]
		set f_db_base 	[string range $f_db_file 0 end-4]

		## If file doesn't already exist, blast in one direction
		if {[file exists "$out_sub/$f_query_base\&$f_db_base\.tsv"] == 0} {
			catch {exec blastp -query $proteome -db "$grp_db_dir/$f_db_base" -out "$out_sub/$f_query_base\&$f_db_base\.tsv" -evalue 1e-10 -outfmt 6 -max_target_seqs 1 -max_hsps 1 -seg "yes" -soft_masking true -num_threads 1}
			set blast_run true			
		}

		## Get the proteome name from the previous db_name
		set r_query_file	"$f_db_base\_proteome.fasta"
		set r_db_base		$f_query_base

		## If file doesn't already exist, then blast in the other direction
		if {[file exists "$out_sub/$f_db_base\&$f_query_base\.tsv"] == 0} {
			catch {exec blastp -query "$grp_in_dir/$r_query_file" -db "$all_db_dir/$r_db_base" -out "$out_sub/$f_db_base\&$f_query_base\.tsv" -evalue 1e-10 -outfmt 6 -max_target_seqs 1 -max_hsps 1 -seg "yes" -soft_masking true -num_threads 1}
			set blast_run true
		}
			
	}
}

if {$blast_run == true} {
	puts stdout		"BLAST was run"
} else {
	puts stdout		"BLAST was not run"
}