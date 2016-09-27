source ~/Dropbox/Scripts/General_utils.tcl
source ~/Dropbox/Scripts/Procs/gc_content.tcl

package require sqlite3

cd /Users/aesin/Desktop/Geo_analysis/Geo_ortholog_nucl
sqlite3 db1 Geo_v_all_orth_location_db

cd /Users/aesin/Desktop/Geo_analysis/Geo_omes/Genome_alignment/Geobac_genomes_processed

set genomes [glob *fna]
set output_table {Genome\tGC_content}

foreach genome $genomes {

	set ncbi_id [string range $genome 0 [string first "\_genomic" $genome]-1]
	
	set binomial [db1 eval {SELECT binomial from t1 WHERE ncbi_id = $ncbi_id LIMIT 1}]
	puts $binomial

	set all_contigs [string trim [openfile $genome]]
	set name [string trim [string range $all_contigs 0 [string first \n $all_contigs]]]

	set gc_cont [genome_gc_cont_all $genome]
	lappend output_table "$binomial\t$gc_cont"

}

set out [open /Users/aesin/Desktop/Geo_analysis/Geo_omes/Geobac_gc_content/Geo_gc_table.tsv w]
puts $out [join $output_table \n]
close $out