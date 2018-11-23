#!/usr/local/bin/tclsh

source ~/Documents/Scripts/General_utils.tcl

set master_dir		/Users/aesin/Desktop/HTa_tonio

set clean_prot_dir	$master_dir/Proteomes/Proteome_clean
set blast_dir		$master_dir/Blast
set blast_in_dir	$blast_dir/HTa_sequence
set blast_db_dir	$blast_dir/Blast_db
set blast_out_dir	$blast_dir/Blast_result

set blast_input_file	$blast_in_dir/HTa_P02345.fasta
set blast_genomes		[glob $clean_prot_dir/*.fasta]

file mkdir $blast_db_dir; file mkdir $blast_out_dir

# Make blast database
foreach genome $blast_genomes {
	set core_name [string range [file tail $genome] 0 end-15]
	if {[file exists $blast_db_dir/$core_name\.phr] != 1} {
		exec makeblastdb -in $genome -dbtype prot -out $blast_db_dir/$core_name
	}
}

# Run the blast
foreach target_genome $blast_genomes {

	set core_name [string range [file tail $target_genome] 0 end-15]

	catch {exec blastp -query $blast_input_file -db $blast_db_dir/$core_name -out $blast_out_dir/$core_name\.tsv -evalue 1e-10 -outfmt 6 -seg "yes" -soft_masking true -num_threads 10}
}
