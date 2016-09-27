#!/usr/local/bin/tclsh

source ~/Dropbox/Scripts/General_utils.tcl
source ~/Dropbox/Scripts/Procs/gc_content.tcl
package require sqlite3

###########################################################################

set penalty 5

###########################################################################
proc geo_origin_positions {} {
	puts "Getting a table of the dnaA / dnaN gene positions to determine genome origin ..."

	set origin_genes_l [split [string trim [openfile /Users/aesin/Desktop/Geo_analysis/HGT_position/Geo_origin_genes.txt]] \n]

	global ori_rep_tbl {}
	global species_l {}

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

	puts "DONE: Table of origin positions"
	return
}

proc geo_assembled_only {} {

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
	puts "DONE: Identified which species are fully assembled genomes"
	return $assembled_species_l
}


proc gc_content_for_segment {dna_sequence} {

	set seq_length [string length $dna_sequence]

	set nG [regsub -all "G" $dna_sequence {} dna_block_mod]
	set nC [regsub -all "C" $dna_block_mod {} dna_block_mod]
	set nN [regsub -all "N" $dna_block_mod {} dna_block_mod]

	set gc_content [tcl::mathfunc::roundto [expr double(($nG) + $nC + [expr double($nN) / 2]) / $seq_length ] 4]

	return $gc_content
}

proc event_position_process {event gapless_dna} {
	set split_by_col [split $event " "]
	global relative_pos

	set start	[lindex $split_by_col 5]
	set end		[lindex $split_by_col 6]
	set strand	[lindex $split_by_col 7]
	set relative_pos [lindex $split_by_col end]
	#puts $strand

	## The genbank positions are 1-indexed, while tcl is 0-indexed ##
	set gene [string range $gapless_dna $start-1 $end-1]
	return $gene
}

###########################################################################

set direct /users/aesin/Desktop/Geo_analysis/Geo_omes/Assembled_geo_no_plasmids

set output_dir /users/aesin/desktop/Geo_analysis/HGT_position/GC_content_per_gene/T$penalty
file mkdir $output_dir/Vert
file mkdir $output_dir/Long

## Open the ortholog database ##
cd /Users/aesin/Desktop/Geo_analysis/Geo_ortholog_nucl
sqlite3 db1 Geo_v_all_orth_location_NEW_db

###############################################################
## Table of species with their dnaA and polB gene locus tags ##
###############################################################
geo_origin_positions

#####################################
## Which genomes are contig-based? ##
#####################################
set assembled_species_l [geo_assembled_only]

##################################################
## Total genome lengths (adjusted for plasmids) ##
##################################################

set genome_length_l [split [string trim [openfile  /Users/aesin/Desktop/Geo_analysis/Geo_omes/Genome_lengths.tsv]] \n]
puts "DONE: Table of genome lengths\n"

####################################################################################################################
## For each assembled genome, get a list of vertical and horizontal GC contents corresponding with their position ##
####################################################################################################################

foreach species $assembled_species_l {

	## Get ncbi_id ##
	set ncbi_id [db1 eval {select ncbi_id from t1 where binomial = $species LIMIT 1}]
	puts $ncbi_id

	## Open the corresponding fasta file (indexed at 0, same as genbank positions). Make gapless ##
	set genome_fasta	[string trim [openfile $direct/Assembled_no_plasmid_fastas/$ncbi_id\_genomic.fasta]]
	set headerless_fasta [string trim [string range $genome_fasta [string first \n $genome_fasta] end]]
	regsub -all "\n" [string toupper $headerless_fasta] {} no_gaps_dna

	## Open the file containing the positions of all the vertical (non-plasmid) genes created in Per_penalty_position.tcl (HGT Position folder) ##
	set vert_events		[split [string trim [openfile $direct/../../HGT_position/For_circular/T$penalty/Per_species_vert/Full_entries/$species\_vert.tsv]] \n]

	set $species\_vert_out {}
	foreach event $vert_events {
		set gene [event_position_process $event $no_gaps_dna]
		set gene_gc_cont [gc_content_for_segment [string trim $gene]]
		
		lappend $species\_vert_out "$relative_pos\t$gene_gc_cont"
	}

	## Open the file containing the positions of all the long HGT (non-plasmid) genes created in Per_penalty_position.tcl (HGT Position folder) ##
	set long_events		[split [string trim [openfile $direct/../../HGT_position/For_circular/T$penalty/Per_species_long/Full_entries/$species\_HGT.tsv]] \n]

	set $species\_long_out {}
	foreach event $long_events {
		set gene [event_position_process $event $no_gaps_dna]
		set gene_gc_cont [gc_content_for_segment [string trim $gene]]
		
		lappend $species\_long_out "$relative_pos\t$gene_gc_cont"
	}

	set vert_out_tbl [join [expr $$species\_vert_out] \n]
	set long_out_tbl [join [expr $$species\_long_out] \n]

	set out [open $output_dir/Vert/$species\_per_vert_gene.tsv w]
	puts $out $vert_out_tbl
	close $out

	set out [open $output_dir/Long/$species\_per_long_gene.tsv w]
	puts $out $long_out_tbl
	close $out
}

db1 close


#db1 close