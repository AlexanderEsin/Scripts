## Blast unaligned contigs to the plasmid sequence database and remove any unaligned contigs that have significant hits to entries in the plasmid database ##

source ~/Dropbox/Scripts/General_utils.tcl

set direct /users/aesin/desktop/Geo_analysis/Geo_omes/Genome_processing/Unaligned_contigs

cd $direct/Fasta
set input_fasta_files [glob *fna]

cd $direct

foreach fasta_file $input_fasta_files {
	set output_name "Blast_out_[string range $fasta_file 0 end-4]\.tsv"
	catch {exec blastn -query ./Fasta/$fasta_file -db /users/aesin/desktop/Blast/Plasmid_seq/DB/plasmid_seq -out ./Blast_v_plasmid_out/$output_name -evalue 1e-50 -outfmt 6 -max_target_seqs 1 -num_threads 4}
}

### Remove any unaligned contigs that had hits to the plasmid database in the blast above. The stringent e-val of -50 should be sufficient to still include weak hits. Although it cannot be excluded that the contigs are in fact genomic and originate from plasmid incorporation into the chromosome, this is less likely - and difficult to determine with certainty ###

puts "Removing contigs that had significant BLAST hits to the plasmid nucleotide database"
set direct /users/aesin/desktop/Geo_analysis/Geo_omes/Genome_processing/

cd $direct/Geobac_genomes_contigs_processing/Geobac_genomes_contig_full_order
set full_ordered_files [glob *fna]
puts "Genome files: [llength $full_ordered_files]"

cd $direct

foreach full_file $full_ordered_files {
	puts "\nWorking on [string range $full_file 0 end-4]"
	set genome [openfile $direct/Geobac_genomes_contigs_processing/Geobac_genomes_contig_full_order/$full_file]
	set all_contigs [split_genes $genome]
	puts "Total number of contigs: [llength $all_contigs]"

	## Get the list of contig with blast hits ##
	set blast_out_file_name "Blast_out_[string range $full_file 0 [string first "genomic" $full_file]-2]\_unaligned_contigs.tsv"
	#puts $blast_out_file_name
	set blast_results [split [openfile $direct/Unaligned_contigs/Blast_v_plasmid_out/$blast_out_file_name] \n]
	puts "Total number of blast hits (may be redundant): [llength $blast_results]"
	
	set contigs_removed {}
	set num_contigs_removed 0
	foreach blast_hit $blast_results {
		set contig_name [lindex [split $blast_hit \t] 0]
		if {[lsearch $contigs_removed $contig_name] == -1} {
			set contig_seq_to_remove [lsearch -glob -inline $all_contigs *$contig_name*]
			set all_contigs [lremove $all_contigs $contig_seq_to_remove]
			incr num_contigs_removed
			lappend contigs_removed $contig_name
		}
	}
	puts "Number of contigs removed: $num_contigs_removed"
	puts "New number of contigs: [llength $all_contigs]"

	set out [open $direct/Geobac_genomes_contigs_processing/Geobac_genomes_contig_trimmed_order/[string range $full_file 0 end-4]_trimmed.fna w]
	puts $out [join $all_contigs ""]
	close $out
}