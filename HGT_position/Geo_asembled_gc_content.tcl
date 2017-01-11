#!/usr/local/bin/tclsh

source ~/Dropbox/Scripts/General_utils.tcl
source ~/Dropbox/Scripts/Procs/gc_content.tcl
package require sqlite3

###############
set window 1000
###############

proc gc_content_by_window {dna_sequence {window 500}} {
	set gc_cont_list {}
	set sequence_length [string length $dna_sequence]

	progress_init $sequence_length
	set counter 0

	set dna_list [split $dna_sequence {}]
	while {$counter < $sequence_length} {
		set segment_l [lrange $dna_list $counter [expr ($counter + $window - 1)]]
		set segment_length [llength $segment_l]
		set segment_str [join $segment_l {}]

		set nG [regsub -all "G" $segment_str {} dna_block_mod]
		set nC [regsub -all "C" $dna_block_mod {} dna_block_mod]
		set nN [regsub -all "N" $dna_block_mod {} dna_block_mod]

		set block_gc_content [tcl::mathfunc::roundto [expr double(($nG) + $nC + [expr double($nN) / 2]) / $segment_length ] 4]
		lappend gc_cont_list $block_gc_content

		set counter [expr $counter + $window]
		progress_tick $counter
	}
	return $gc_cont_list
}

proc reverse_compl {sequence} {

	if {[regexp -all " " $sequence] > 0 || [regexp -all "\n" $sequence] > 0 || [regexp -all "\t" $sequence] > 0} {
	   	puts "There are gaps in the sequence - exiting"
	   	return
	}
	
	set split_by_nucl [split [string toupper $sequence] {}]
	set join_w_space " [string trim [join $split_by_nucl " "]]"
	set counter [regsub -all " A" $join_w_space {T} join_mod]
	set counter [expr $counter + [regsub -all " T" $join_mod {A} join_mod]]
	set counter [expr $counter + [regsub -all " G" $join_mod {C} join_mod]]
	set counter [expr $counter + [regsub -all " C" $join_mod {G} join_mod]]

	set leftover_bases [regexp -all " " $join_mod]
	if {$leftover_bases > 0} {
		puts "There are $leftover_bases non-canonical bases in this sequence. Changing them to N."
		regsub -all {\s{1}[A-Z]{1}} $join_mod {N} join_mod
	}

	set reversed_seq [string_reverse $join_mod]
	return $reversed_seq
	
}

###########################################################################
set direct /users/aesin/Desktop/Geo_analysis/Geo_omes/Assembled_geo_no_plasmids

set output_dir /users/aesin/desktop/Geo_analysis/HGT_position/GC_content/Window_$window
file mkdir $output_dir

## Open the ortholog database ##
cd /Users/aesin/Desktop/Geo_analysis/Geo_ortholog_nucl
sqlite3 db1 Geo_v_all_orth_location_NEW_db

###############################################################
## Table of species with their dnaA and polB gene locus tags ##
###############################################################

puts "Getting a table of the dnaA / dnaN gene positions to determine genome origin ..."

set origin_genes_l [split [string trim [openfile /Users/aesin/Desktop/Geo_analysis/HGT_position/Geo_origin_genes.txt]] \n]

set ori_rep_tbl {}
set species_l {}

foreach entry $origin_genes_l {
	set entry_l [split $entry \t]
	set species [lindex $entry_l 0]

	set dnaA_gene [lindex $entry_l 1]
	set dnaA_start_end [db1 eval {select gene_start, gene_end, strand from t1 where locus_tag = $dnaA_gene}]

	## If there's no entry for dnaA, resort to polB (dnaN) ##
	if {[string length $dnaA_start_end] == 0} {
		set polB_gene [lindex $entry_l 2]
		set polB_start_end [db1 eval {select gene_start, gene_end, strand from t1 where locus_tag = $polB_gene}]

		set polb_position_zero	[lindex $polB_start_end 0]
		set polb_position_last	[lindex $polB_start_end 1]
		set polb_strand			[lindex $polB_start_end 2]
		lappend ori_rep_tbl "$species\t$polb_position_zero\t$polb_position_last\t$polb_strand"
	} else {
		set dnaa_position_zero	[lindex $dnaA_start_end 0]
		set dnaa_position_last	[lindex $dnaA_start_end 1]
		set dnaa_strand			[lindex $dnaA_start_end 2]
		lappend ori_rep_tbl "$species\t$dnaa_position_zero\t$dnaa_position_last\t$dnaa_strand"
	}

	lappend species_l $species
}

set out [open /Users/aesin/Desktop/Geo_analysis/HGT_position/Geo_origin_positions.tsv w]
puts $out [join $ori_rep_tbl \n]
close $out
puts [join $ori_rep_tbl \n]
puts "DONE: Table of origin positions\n"


##################################################
## Total genome lengths (adjusted for plasmids) ##
##################################################

set genome_length_l [split [string trim [openfile  /Users/aesin/Desktop/Geo_analysis/Geo_omes/Genome_lengths.tsv]] \n]
puts "DONE: Table of genome lengths\n"

#####################################
## Which genomes are contig-based? ##
#####################################

if {[file exists /users/aesin/desktop/Geo_analysis/HGT_position/assembled_geo.tsv] == 1} {
	set assembled_species_l [split [string trim [openfile /users/aesin/desktop/Geo_analysis/HGT_position/assembled_geo.tsv]] \n]
} else {
	cd /Users/aesin/Desktop/Geo_analysis/Geo_omes/Genome_processing/Geobac_genomes_assembled_genbank
	set non_conting_genomes [glob *gbk]
	regsub -all {_genomic.gbk} $non_conting_genomes {} non_contig_ncbi_ids

	set assembled_species_l {}
	foreach assembled_ncbi_id $non_contig_ncbi_ids {
		set assembled_species [db1 eval {select binomial from t1 where ncbi_id = $assembled_ncbi_id LIMIT 1}]
		lappend assembled_species_l $assembled_species
	}

	set out [open /users/aesin/desktop/Geo_analysis/HGT_position/assembled_geo.tsv w]
	puts $out [join $assembled_species_l \n]
	close $out
}
puts "DONE: Identified which species are fully assembled genomes\n"

##################################################################
## Create GBK files which are lacking plasmids (assembled only) ##
##################################################################

set partition_pat {\n+(//)+\n+}
cd $direct/All_geobac_gbk

set geo_gbk_files [glob *gbk]
puts "Working with [llength $geo_gbk_files] Geobac .gbk files\n"

foreach gbk_file $geo_gbk_files {

	# Get the binomial name #
	set ncbi [string range $gbk_file 0 [string last "_" $gbk_file]-1]
	set binomial [db1 eval {select binomial from t1 where ncbi_id = $ncbi LIMIT 1}]

	# Only continue if the genome is assembled #
	if {[lsearch $assembled_species_l $binomial] == -1} {
		continue
	}

	set gbk_data [string trim [openfile $gbk_file]]
	
	regsub -all $partition_pat $gbk_data "££" genome_data
	set partitions [split $genome_data {££}]

	set chromosome_only {}

	foreach partition $partitions {

		if {[string length $partition] == 0} {
			continue
		}
		set DEFINITION_line [regexp -line -inline {(?:DEFINITION)+.+$} $partition]
		if {[regexp {plasmid} $DEFINITION_line] == 0} {
			lappend chromosome_only $partition
		}
	}

	set new_gbk [join $chromosome_only "\n//\n"]
	set out [open $direct/Assembled_no_plasmid_gbks/$gbk_file w]
	puts $out $new_gbk
	close $out

	puts "Genome: $binomial\t$ncbi\tPost-processing partitions: [llength $chromosome_only]"
}

############################################################################
## Convert these non-plasmid GBK files into fasta files with DNA sequence ##
############################################################################

cd $direct/Assembled_no_plasmid_gbks
set gbk_files [glob *gbk]

foreach gbk_file $gbk_files {
	set out_fasta_name "[string range $gbk_file 0 end-4]\.fasta"
	catch {exec seqret -snucleotide1 Yes $gbk_file $direct/Assembled_no_plasmid_fastas/$out_fasta_name -osformat fasta}
}

#########################################################
## Rearrange the fasta genomes to start at the origin  ##
#########################################################
cd $direct/Assembled_no_plasmid_fastas
set fasta_files [glob *fasta]
puts stdout "\nRearranging genome sequence to put origin at position 0 and correct for orientation...\n"

foreach fasta_file $fasta_files {
	set ncbi [string range $fasta_file 0 [string last "_" $fasta_file]-1]
	set binomial [db1 eval {select binomial from t1 where ncbi_id = $ncbi LIMIT 1}]

	set relevant_genome_length [lindex [split [lsearch -inline -glob $genome_length_l $binomial\t*] \t] 1]

	set species_start [lsearch -glob -inline $ori_rep_tbl $binomial\t*]

	set genome [string trim [openfile $fasta_file]]
	set header [string trim [string range $genome 0 [string first \n $genome]]]

	regsub -all {\n} [string trim [string range $genome [string first \n $genome] end]] {} dna_no_gaps
	puts "Species: $binomial\t\tGenome length: [string length $dna_no_gaps]"

	set ori_strand [lindex [split $species_start \t] 3]
	if {$ori_strand eq "+"} {
		set ori_start [lindex [split $species_start \t] 1]
		set downstream [string range $dna_no_gaps $ori_start-1 end]
		set upstream [string range $dna_no_gaps 0 $ori_start-2]

		set recombine [append downstream $upstream]

	} else {
		set rev_comp [reverse_compl $dna_no_gaps]
		set ori_start [expr [string length $rev_comp] - [lindex [split $species_start \t] 2]]

		set downstream [string range $rev_comp $ori_start end]
		set upstream [string range $rev_comp 0 $ori_start-1]

		set recombine [append downstream $upstream]
	}

	puts "Rearrangment complete.\t\tLength: [string length $recombine]"
	set out [open $direct/Assembled_no_plasmid_rearranged_fastas/$ncbi\_rearranged.fasta w]
	puts $out [string toupper $recombine]
	close $out

}

##########################################################
## Process the files and perform GC content calculation ##
##########################################################

cd $direct/Assembled_no_plasmid_rearranged_fastas
set fasta_files [glob *fasta]
puts "\n\nCalculating GC content based on a window size of: $window"

foreach fasta_file $fasta_files {
	set ncbi [string range $fasta_file 0 [string last "_" $fasta_file]-1]
	set binomial [db1 eval {select binomial from t1 where ncbi_id = $ncbi LIMIT 1}]

	set genome [string trim [openfile $fasta_file]]
	set genome_length [string length $genome]
	
	puts "Working on $binomial ..."
	set gc_content_l [gc_content_by_window $genome $window]

	## Prepare the output ##

	set half_window [expr $window / 2]
	set coordinate $half_window

	set remainder [expr $genome_length % $window]
	#puts $remainder

	set output_coord_gc_tbl {}
	set i 0

	foreach gc_measure $gc_content_l {
		incr i
		set fraction_coord [tcl::mathfunc::roundto [expr double($coordinate) / double($genome_length)] 5]

		set coord_gc_entry "$fraction_coord\t$gc_measure"
		lappend output_coord_gc_tbl $coord_gc_entry
		## If this is the penultimate loop, we don't want to add a whole window to the coordinate - just the remainder of the genome ##
		if {$i == [expr [llength $gc_content_l] - 1]} {
			set coordinate [expr $coordinate + $remainder]
		} else {
			set coordinate [expr $coordinate + $window]
		}		
	}

	set out [open $output_dir/$binomial\_gc_cont.tsv w]
	puts $out [join $output_coord_gc_tbl \n]
	close $out
}

puts "\nDONE"


db1 close

puts stdout $old_orient_dnaa_gene




















