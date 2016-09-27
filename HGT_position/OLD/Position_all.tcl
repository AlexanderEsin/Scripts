#!/usr/local/bin/tclsh

source ~/Dropbox/Scripts/General_utils.tcl
package require sqlite3

###########################################################################
## If linear representation, then linear == TRUE, if circular linear == FALSE

set linear FALSE

###########################################################################
if {$linear == TRUE} {
	set output_dir /users/aesin/desktop/Geo_analysis/HGT_position/For_linear/Per_species_all
} else {
	set output_dir /users/aesin/desktop/Geo_analysis/HGT_position/For_circular/Per_species_all
}

file mkdir $output_dir/Full_entries
file mkdir $output_dir/Start
file mkdir $output_dir/Relative_start

## Open the ortholog database ##
cd /Users/aesin/Desktop/Geo_analysis/Geo_ortholog_nucl
sqlite3 db1 Geo_v_all_orth_location_NEW_db

###########################################################################
## Table of species with their dnaA and polB gene locus tags ##
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

puts "DONE: Table of origin positions"

###########################################################################
## Total genome lengths (adjusted for plasmids) ##
set genome_length_l [split [string trim [openfile  /Users/aesin/Desktop/Geo_analysis/Geo_omes/Genome_lengths.tsv]] \n]
puts "DONE: Table of genome lengths"

###########################################################################
## Which genomes are contig-based? ##
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
###########################################################################

set plasmid_gene 0

set plasmid_entries {}
set species_list {}

set group_counter 0
set total_genes_counter 0
set total_genes_added_to_table 0
set genes_added_assembled 0

###########################################################################
set geo_groups [db1 eval {select DISTINCT group_number from t1 where binomial glob "Geobacillus*"}]


foreach group $geo_groups {
	
	set geo_ids [db1 eval {SELECT prot_id from t1 where group_number = $group AND binomial glob "Geobacillus*";}]
	
	#########################################################################
	## For each transfer event establish the recipient Geobacillus species ##
	#########################################################################

	foreach geo_id $geo_ids {

		incr total_genes_counter
		
		set entry_l [split [db1 eval {SELECT * from t1 where prot_id = $geo_id}] " "]
		if {[llength $entry_l] == 0} {
			incr failed_to_map
			continue
		}

		set species [lindex $entry_l 9]
		# Only continue if the genome is assembled #
		if {[lsearch $assembled_species_l $species] == -1} {
			continue
		}

		################################################################################
		## Check the protein ID against the plasmid TAG list for the relevent species ##
		################################################################################

		## Open Plasmid tag file (if one exists for this species) ##
		set plasmid_tag_file /users/aesin/Desktop/Geo_analysis/Geo_omes/Plasmid_tags/$species\_plasmid_locus_tags.tsv
		if {[file exists $plasmid_tag_file] == 1} {
			set tag_list [split [string trim [openfile $plasmid_tag_file]] \n]
			set transfer_locus [lindex $entry_l 8]

			## If this gene DOES map to the plasmid, do NOT add it to the total output list ##
			if {[lsearch $tag_list $transfer_locus] != -1} {
				incr plasmid_gene
				lappend plasmid_entries $entry_l
				continue
			}
		}

		######################################################################################
		## Put the information (either full or just start position) into per-species tables ##
		######################################################################################

		## Full table ##
		if {[info exists $species\_tbl] == 0} {
			set $species\_tbl {}
			lappend species_list $species
		}

		## The strand will determine whether the true gene start position is the start or the end position as in the table. I.e. for complement strand, the "start" in the table is in fact the end, and vice versa ##
		## Start positions only ##
		set strand [lindex $entry_l 7]
		if {$strand eq "+"} {
			set start_pos [lindex $entry_l 5]
		} else {
			set start_pos [lindex $entry_l 6]
		}
				
		if {[info exists $species\_start_tbl] == 0} {
			set $species\_start_tbl {}
		} 
		lappend $species\_start_tbl $start_pos

		###########################################################################
		## Make the distance to the origin a relative value of the genome length ##
		###########################################################################

		## Relative start positions ##
		if {[regexp "Geobacillus" $species] == 1} {
			set relevant_genome_length [lindex [split [lsearch -inline -glob $genome_length_l $species\t*] \t] 1]

			## Get the origin position for this genome ##
			set species_start [lsearch -glob -inline $ori_rep_tbl $species\t*]

			## If the origin lies on the complement strand, we take the "end" as the start. We also change the upstream/downstream orientation. ##
			set ori_strand [lindex [split $species_start \t] 3]
			if {$ori_strand eq "+"} {
				set ori_start [lindex [split $species_start \t] 1]

				# Get the relative position #
				set relative_pos [expr $start_pos - $ori_start]
				# If it's negative then it lies upstream of the origin and we append it to the "end" #
				if {$relative_pos < 0} {
					set relative_pos [expr $relative_pos + $relevant_genome_length]
				}

			} else {
				set ori_start [lindex [split $species_start \t] 2]

				# Get the relative position #
				set relative_pos [expr $ori_start - $start_pos]
				# If it's negative then it lies upstream of the origin and we append it to the "end" #
				if {$relative_pos < 0} {
					set relative_pos [expr $relative_pos + $relevant_genome_length]
				}

			}
			
			## Get the fractional position ##
			set fraction_pos [tcl::mathfunc::roundto [expr double($relative_pos) / double($relevant_genome_length)] 5]
			## Append the fractional position to the full entry table ##
			lappend $species\_tbl [join [concat $entry_l $fraction_pos] " "]

			## Add it to the per-species table ##
			if {[info exists $species\_rel_start_tbl] == 0} {
				set $species\_rel_start_tbl {}
			} 
			lappend $species\_rel_start_tbl $fraction_pos

			## If the position is in the first third or so of the genome, add 1 to have it rollover so that we can see the correct density across the origin ##
			if {$linear == TRUE && $fraction_pos < 0.3} {
				set rollover_fraction_pos [expr $fraction_pos + 1]
				lappend $species\_rel_start_tbl $rollover_fraction_pos	
			}
		}

	}

	incr group_counter
	puts "Groups processed: $group_counter / [llength $geo_groups]"
	
}

set combined_rel_pos {}
foreach species $species_list {
	set full_out_tbl [join [expr $$species\_tbl] \n]
	set start_out_tbl [join [expr $$species\_start_tbl] \n]
	

	set out [open $output_dir/Full_entries/$species\_all.tsv w]
	puts $out $full_out_tbl
	close $out

	set out [open $output_dir/Start/$species\_all.tsv w]
	puts $out $start_out_tbl
	close $out

	if {[regexp "Geobacillus" $species] == 1} {
		set rel_start_out_tbl [join [expr $$species\_rel_start_tbl] \n]
		set out [open $output_dir/Relative_start/$species\_all.tsv w]
		puts $out $rel_start_out_tbl
		close $out

		if {[lsearch $assembled_species_l $species] != -1} {
			puts "\nAdding species $species to table because it's fully assembled"
			lappend combined_rel_pos $rel_start_out_tbl
			set genes_added_assembled [expr $genes_added_assembled + [llength $rel_start_out_tbl]]
		} else {
			puts "$species not added to the table as it's contig based. It has [llength $rel_start_out_tbl] position data."
		}
	}	
}

set out [open $output_dir/combined_rel_pos_assembled.tsv w]
puts $out [join $combined_rel_pos \n]
close $out

db1 close

puts "Total num genes processed: $total_genes_counter\n"
puts "Total num genes added to table: $total_genes_added_to_table\n"
puts "Genes added from fully assembled species: $genes_added_assembled"
puts "Genes mapping to plasmid: $plasmid_gene"
#puts "Plasmid entries:\n[join $plasmid_entries \n]"