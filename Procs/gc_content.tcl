source ~/Dropbox/Scripts/General_utils.tcl


proc reform_fna_fasta {fasta_file} {
	global header_db_entry
	global no_gaps_fasta
	global element_id_list

	set output {}
	set header_db_entry {}
	set element_id_list {}

	set unformatted_sequence [openfile $fasta_file]
	set genome_elements [split_genes $unformatted_sequence]

	foreach element $genome_elements {
		set element [string trim $element]
		set header [string trim [string range $element 0 [string first \n $element]]]
		set element_id [string trim [string range $header 1 [string first " " $header]]]

		lappend element_id_list $element_id
		lappend header_db_entry "$element_id\t[string trim [string range $header [string first " " $header] end]]"

		regsub -all "\n" [string range $element [string first \n $element] end] {} element_seq_no_gaps
		lappend output ">$element_id\n$element_seq_no_gaps"
	}

	set element_id_list [join $element_id_list \n]
	set header_db_entry [join $header_db_entry \n]
	set no_gaps_fasta [join $output \n]
}

### Updated much faster version in the Geo_assembled_gc_content.tcl script ###
# # Reverse complement given nucleotide sequence #
# proc reverse_compl {nuc_sequence} {
#     global rc_sequence

#     set sequence [string trim $nuc_sequence]

#     if {[regexp -all " " $sequence] > 0 || [regexp -all "\n" $sequence] > 0} {
#     	puts "There are gaps in the sequence - exiting"
#         return
#     } else {
#         set split_by_nucl [split $sequence {}]
#         progress_init [llength $split_by_nucl]
#         set counter 0
#         set rc_sequence {}
#         foreach nucl $split_by_nucl {
#             if {$nucl eq "A"} {
#                 set rc_sequence [lreplace $rc_sequence 0 -1 "T"]
#             } elseif {$nucl eq "T"} {
#                 set rc_sequence [lreplace $rc_sequence 0 -1 "A"]
#             } elseif {$nucl eq "G"} {
#                 set rc_sequence [lreplace $rc_sequence 0 -1 "C"]
#             } elseif {$nucl eq "C"} {
#                 set rc_sequence [lreplace $rc_sequence 0 -1 "G"]
#             } elseif {$nucl eq "N" || $nucl eq "X"} {
#                 set rc_sequence [lreplace $rc_sequence 0 -1 "N"]
#             ## If the nomenclature is not the standard bases + N ... ##
#             } else {
#             	# If puRine convert to pYrimidine #
#             	if {$nucl eq "R"} {
#             		set rc_sequence [lreplace $rc_sequence 0 -1 "Y"]
#             	# If pYrimidine convert to puRine #
#             	} elseif {$nucl eq "Y"} {
#             		set rc_sequence [lreplace $rc_sequence 0 -1 "R"]
#             	# If any other character -> conver to N #
#             	} else {
#             		set rc_sequence [lreplace $rc_sequence 0 -1 "N"]
#             	}      	
#             }

#             progress_tick $counter
#             incr $counter
#         }
#         set rc_sequence [join $rc_sequence {}]
#         return $rc_sequence
#     }
# }

## Calculate the GC content of a sequence, given as a string (can be fast formatted but without header) ###
proc gene_gc_cont {gene_string} {
	global gene_length
	global gc_content

	regsub -all "\n" [string trim $gene_string] {} gene_string
	set gene_length [string length $gene_string]

	# If there are any non-standard bases, report their number #
	set num_non_standard [regexp -all {[^AGTC]} $gene_string]
	if {$num_non_standard != 0} {
		set removed [regsub -all "N" $gene_string {} gene_string]
		puts "Sequence length: $gene_length\tNon-standard bases: $num_non_standard\tRemoved Ns: $removed"
	}

	# Count the G/C number and also S bases (Strong) #
	set nG [regsub -all "G" $gene_string {} gene_string_mod]
	set nC [regsub -all "C" $gene_string_mod {} gene_string_mod]
	set nStr [regsub -all "S" $gene_string_mod {} gene_string_mod]

	# GC content is Sum(G + C + Str)/Nucleotide_length #
	set gc_content [expr roundto((double($nG) + $nC + $nStr) / $gene_length)]
	return $gc_content
}

## NEEDS TO BE REWORKED ALONG THE LUNES OF gene_gc_cont ##
proc genome_gc_cont_all {fasta_file} {
	global genome_length
	global gc_content

	set all_contigs [openfile $fasta_file]
	set genome_elements [regsub -all -line {>.+} [string trim $all_contigs] {} dna_seq]
	regsub -all "\n" $dna_seq {} gapless_dna_seq

	set genome_length [string length $gapless_dna_seq]
	set gapless_dna_seq [string toupper $gapless_dna_seq]

	set num_non_standard [regexp -all {[^AGTC]} $gapless_dna_seq]
	if {$num_non_standard != 0} {
		set removed [regsub -all "N" $gapless_dna_seq {} gapless_dna_seq]
		puts "Sequence length: $genome_length\tNon-standard bases: $num_non_standard\tRemoved Ns: $removed"
	}

	set nG [regsub -all "G" $gapless_dna_seq {} gapless_dna_seq_mod]
	set nC [regsub -all "C" $gapless_dna_seq_mod {} gapless_dna_seq_mod]
	set nStr [regsub -all "S" $gapless_dna_seq_mod {} gapless_dna_seq_mod]

	set gc_content [expr roundto((double($nG) + $nC + $nStr) / $genome_length)]
	return $gc_content
}

proc new_id_from_db {old_id} {
	global new_id
	set new_id ""
	db1 eval {SELECT New_id FROM t1 WHERE Old_id = $old_id} {
	    set new_id "$New_id"
	}
	return $new_id
}

proc features_from_tbl {feat_tbl_hit} {
	global accession
	global product_name
	global strand
	global coordinate_1
	global coordinate_2
	
	set features [split [lindex $feat_tbl_hit 0] \t]
	set accession [lindex $features 6]
	set coordinate_1 [lindex $features 7]
	set coordinate_2 [lindex $features 8]
	set strand [lindex $features 9]
	set product_name [lindex $features 11]
}

proc get_locus_tag_from_protein_header {protein_hit} {
	global locus_tag

	set header_features [split [string trim [string range $protein_hit 0 [string first \n $protein_hit]]] \t]
	set locus_tag [lindex $header_features 2]
	return $locus_tag
}

proc features_from_protein_header {protein_hit} {
	global accession
	global strand
	global gene_start
	global gene_end

	set header_features [split [string trim [string range $protein_hit 0 [string first \n $protein_hit]]] \t]
	# Process the accession name #
	set accession [lindex $header_features 3]
	if {[regexp {:} $accession] != 0} {
		set accession [lindex [split $accession \:] 0]
	}
	# Process the location data #
	set location [lindex $header_features 4]
	set location_split [split $location " "]
	if {[string trim [lindex $location_split 1]] eq "Forward"} {
		set strand "+"
	} else {
		set strand "-"
	}
	set coordinates [split [lindex $location_split 0] \:]
	set gene_start [lindex $coordinates 0]
	set gene_end [lindex $coordinates 1]
}

proc quick_blast {protein ncbi_id db_path temp_dir} {
	## Prepare the query file with the sequence of interest ##
	set out [open $temp_dir/Temp_query w]; puts $out [string trim $protein]; close $out
	## Set the path to the database ##
	set db_name "$db_path/$ncbi_id\_genomic_db"
	## Execute the blastp and write out to temp directory ##
	catch {exec blastp -query $temp_dir/Temp_query -db $db_name -out $temp_dir/Temp_out -evalue 1e-10 -outfmt 6 -max_target_seqs 1 -max_hsps 1 -num_threads 4}
}

proc quick_process_blast {temp_dir {warnings off}} {
	global new_id_from_blast

	set new_id_from_blast ""
	set results [split [string trim [openfile $temp_dir/Temp_out]] \n]
	if {[llength $results] == 1} {
		set data [split [lindex $results 0] \t]
		set perc_ident [lindex $data 2]
		if {$perc_ident >= 99.00} {
			set old_id [lindex $data 0]
			set new_id_from_blast [lindex $data 1]
			#puts "Found that $old_id MATCHES $new_id exactly ..."
			return $new_id_from_blast
		} else {
			if {$warnings eq "on" || $warnings eq "ON"} {
				puts "Not conclusive --> Percentage identity < 99% \($perc_ident\)"
			}			
			return
		}
	} else {
		if {$warnings eq "on" || $warnings eq "ON"} {
			puts "Not conclusive --> Num. results: [llength $results]"
		}
		return
	}
}
