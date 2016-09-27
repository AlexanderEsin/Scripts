#!/usr/local/bin/tclsh

source ~/Dropbox/Scripts/General_utils.tcl
source ~/Dropbox/Scripts/Procs/gc_content.tcl

package require sqlite3

###########################################################################
set output_folder /Users/aesin/Desktop/Geo_analysis/Geo_omes/Orthocluster_synteny

###########################################################################

## Open the ortholog database ##
cd /Users/aesin/Desktop/Geo_analysis/Geo_ortholog_nucl
sqlite3 db1 Geo_v_all_orth_location_NEW_db

## Get a list of all the groups in the database ##
set groups [lsort -dictionary [db1 eval {SELECT DISTINCT group_number from t1}]]
progress_init [llength $groups]

## Get the table with the total number of Geobacillus sequences in the fasta files of the groups from the stats files ##
set geo_number_table [lrange [split [string trim [openfile /Users/aesin/Desktop/Geo_analysis/Geo_v_all/2.0/Stats/Combined_group_fastas_stats.tsv]] \n] 1 end]

#######################################

## Prepare a list of tables - one for each of the Geobacillus species ##
set geo_species [db1 eval {SELECT DISTINCT binomial from t1 where binomial LIKE 'Geobacillus%'}]
foreach species $geo_species {
	set $species\_table {}
}

## Prepare a list to hold the orthologous family prot_ids ##
set ortholog_id_l {}
## Prepare a list of groups for which we succesfully got positional data ##
set success_group_l {}

#######################################
## Prepare the contig order list for each organism ##
# The all_contigs_l is there to make sure that we are not picking up genes that match removed contigs #
set all_contigs_l {}

file mkdir $output_folder/Contig_order_files

cd /Users/aesin/Desktop/Geo_analysis/Geo_omes/Genome_alignment/Geobac_genomes_processed
set genome_list [glob *fna]
foreach species $geo_species {
	set ncbi_id [db1 eval {SELECT ncbi_id from t1 where binomial = $species LIMIT 1}]
	set genome_file [lsearch -glob -inline $genome_list $ncbi_id*]
	set contigs [split_genes_fast [openfile $genome_file]]

	set contig_order_l [list $species]

	foreach contig $contigs {
		set contig_accession [string range $contig 1 [string first " " $contig]-1]
		lappend contig_order_l $contig_accession
		lappend all_contigs_l $contig
	}

	set out [open $output_folder/Contig_order_files/$species\_contigs.tsv w]
	puts $out [join $contig_order_l \t]
	close $out
}


#######################################


set all_present 0
set not_equal_num_geobac 0
set groups_processed 0

# set groups [lrange $groups 0 0]

foreach group $groups {

	## How many Geobacillus are there in this group in the database ##
	set geo_in_db [db1 eval {SELECT binomial from t1 where group_number = $group AND binomial LIKE 'Geobacillus%'}]
	set num_geo_in_db [llength $geo_in_db]

	## How many Geobacillus are there in this group in our MSA ##
	set num_geo_in_group [lindex [split [lsearch -inline -glob $geo_number_table $group\t*] \t] 2]

	if {$num_geo_in_db != $num_geo_in_group} {
		incr not_equal_num_geobac
		continue
	}

	#######################################

	## Assuming all the Geobacillus sequences in the family are mappable, extract their coordinates ##
	foreach species $geo_in_db {
		set entry [db1 eval {SELECT prot_id, accession, gene_start, gene_end, strand from t1 where group_number = $group AND binomial = $species}]
		regsub -all {\+} $entry {1} entry
		regsub -all {\-} $entry {-1} entry

		## Check for multiple proteins from the same species in a group ##
		set indiv_values [split $entry " "]
		if {[llength $indiv_values] > 5} {
			set processed 0
			while {$processed < [llength $indiv_values]} {
				set unique_entry [join [lrange $indiv_values $processed $processed+4]]
				lappend "$species\_table" $unique_entry
				set processed [expr $processed + 5]
			}
		} else {
			lappend "$species\_table" $entry
		}
	}

	#######################################

	## First round - complete a single set of 19 protein IDs ##
	set prime_id_l {}
	set second_round_tbl {}

	set species_index 0
	foreach species $geo_species {
		set geo_prot_id [db1 eval {SELECT prot_id from t1 where group_number = $group AND binomial = $species}]

		if {[string length $geo_prot_id] > 0} {
			## Are there multiple proteins in this gene family from the same organism (duplication)? ##
			set geo_ids [split $geo_prot_id " "]

			set geo_prot_prime_id [lindex $geo_ids 0]
			lappend prime_id_l $geo_prot_prime_id

			if {[llength $geo_ids] > 1} {
				set orthologs [join [lrange $geo_ids 1 end]]
				set table_entry [join [concat $species_index $orthologs] \t]
				lappend second_round_tbl $table_entry
			}
		} else {
			lappend prime_id_l {}
		}

		incr $species_index
	}

	## Add the primal_id_l to output table ##
	lappend ortholog_id_l [join $prime_id_l \t]

	## Append a line substituting in the other orthologs, if presen ##
	if {[llength $second_round_tbl] > 0} {

		foreach entry $second_round_tbl {
			set line_elements [split $entry \t]
			set index [lindex $line_elements 0]
			set orthologs [split [lindex $line_elements 1]]

			foreach ortholog $orthologs {
				set a_new_line $prime_id_l
				set a_new_line [lreplace $a_new_line $index $index $ortholog]
				lappend ortholog_id_l [join $a_new_line \t]
			}
		}
	}


	#puts [join $ortholog_id_l \n]

	incr groups_processed
	progress_tick $groups_processed

}

file mkdir $output_folder/Genome_coordinate_files

foreach species $geo_species {
	set out [open $output_folder/Genome_coordinate_files/$species\_coordinates.tsv w]
	puts $out [join [expr $$species\_table] \n]
	close $out
}

set out [open $output_folder/Ortholog_id_list.tsv w]
puts $out [join $ortholog_id_l \n]
close $out

# set out [open $output_folder/Groups_extracted.tsv w]
# puts $out [join $success_group_l \n]
# close $out

db1 close