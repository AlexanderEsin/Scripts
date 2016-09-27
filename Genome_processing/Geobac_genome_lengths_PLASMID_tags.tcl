#!/usr/local/bin/tclsh

source ~/Dropbox/Scripts/General_utils.tcl
package require sqlite3

###########################################################################

## Open the ortholog database ##
cd /Users/aesin/Desktop/Geo_analysis/Geo_ortholog_nucl
sqlite3 db1 Geo_v_all_orth_location_NEW_db

## Directories ##
set direct /users/aesin/Desktop/Geo_analysis/Geo_omes
set gbk_dir $direct/Genome_alignment/Geobac_genomes_gbff
set plasmid_tag_out $direct/Plasmid_tags

## Get the input gbk files ##
cd $gbk_dir
set input_gbk_genomes [glob *gbk]

## Split the gbk file into partitions (e.g. either contigs or chromosome/plasmid) ##
set pat {\n+(//)+\n+}

## Set output file ##
set output_tbl {}

foreach genome $input_gbk_genomes {
	set genome_data [string trim [openfile $genome]]
	regsub -all $pat $genome_data "££" genome_data
	set partitions [split $genome_data {££}]

	set total_length 0
	set total_clean_tags {}
	set plasmid_stats_file {}

	## Check the partitions ##
	foreach partition $partitions {
		if {[string length $partition] == 0} {
			continue
		}
		set partition [string trim $partition]

		##################################################################
		## Get the locus line, get length, remove it from the partition ##
		##################################################################
		set LOCUS_line [string range $partition 0 [string first "\n" $partition]]
		# puts $LOCUS_line
		
		regsub -all {\s+} $LOCUS_line "\t" LOCUS_line
	
		set partition [string trim [string range $partition [string first "\n" $partition] end]]

		########################################################################################
		## Get the definition line. Check if plasmid. Only add to total length if NOT plasmid ##
		########################################################################################
		set DEFINITION_line [string range $partition 0 [string first "\n" $partition]]
		if {[regexp -all {plasmid} $DEFINITION_line] == 1} {

			lappend plasmid_stats_file "\n$LOCUS_line"
			lappend plasmid_stats_file $DEFINITION_line

			set locus_tags [regexp -all -inline -line {(?:/locus_tag=)+.+} $partition]
			set clean_partition_tags {}
			
			foreach tag $locus_tags {
				set tag_trim [string range $tag [string first "\"" $tag]+1 [string last "\"" $tag]-1]
				lappend clean_partition_tags $tag_trim
			}

			set clean_partition_tags [lsort -unique $clean_partition_tags]
			lappend plasmid_stats_file [join $clean_partition_tags \n]
			lappend total_clean_tags [join $clean_partition_tags \n]
		} else {
			set total_length [expr $total_length + [lindex [split $LOCUS_line \t] 2]]
		}
	}


	## Get the binomial ##
	set ncbi [string range $genome 0 [string last "_" $genome]-1]
	set binomial [db1 eval {select binomial from t1 where ncbi_id = $ncbi LIMIT 1}]

	## Append the total length of the genome to the global output table ##
	lappend output_tbl "$binomial\t$total_length"

	## Write out the plasmid tags ##
	if {[llength $total_clean_tags] > 0} {
		set out [open $plasmid_tag_out/$binomial\_plasmid_locus_tags.tsv w]
		puts $out [join $total_clean_tags \n]
		close $out

		set out [open $plasmid_tag_out/$binomial\_plasmid_stats.tsv w]
		puts $out [string trim [join $plasmid_stats_file \n]]
		close $out
	}
	
}

## Write out the genome lengths ##
set out [open $direct/Genome_lengths.tsv w]
puts $out [join $output_tbl \n]
close $out


