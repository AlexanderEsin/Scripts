#!/usr/local/bin/tclsh
## Process the output of query fasta vs eggNOG_bact database hmmscan ##
source ~/Dropbox/Scripts/General_utils.tcl
source ~/Dropbox/Scripts/Procs/EggNOG.tcl

###########################################################################

set direct /users/aesin/desktop
set working_directory $direct/Geo_analysis/Geo_ortholog_nucl/Geobac_proteins_no_coords

## Set file names ##
set hmm_out_file All_unmapped_bactNOG_hmmscan_out.txt
set hmm_query_file Geobac_not_mapped_fasta.faa
set parse_process_out_name All_unmapped; # Non_part_hyp_unmapped

## Let's check how many input queries we had ##
set query_input_path $working_directory/Fastas/$hmm_query_file
set input_queries [split_genes [openfile $query_input_path]]
set num_input_queries [llength $input_queries]

## Get a list of protein IDs from the query file ##
set query_ids {}
foreach input $input_queries {
	set id [string range $input 1 [string first " " $input]-1]
	lappend query_ids $id
}

cd $working_directory/Functional_annotation

## Open a log file to which we write some stats ##
set log [open $parse_process_out_name\_stats.txt w]
puts $log "Query fasta input into Hmmscan vs the bactNOG database: $hmm_query_file\n\tNumber of input proteins sequences: $num_input_queries"

## Let's parse the hmmscan output file and identify top eggNOG hits for each query (those that did not find any would not appear in the table) ##
set hmm_out_path $working_directory/Functional_annotation/$hmm_out_file
set output_table [parse_nog_search_output $hmm_out_path]
set num_found_hmm_cluster [llength $output_table]
puts $log "\nHmmscan output file name: $hmm_out_file\n\tNumber of sequences finding hits in bactNOG database: $num_found_hmm_cluster"

## The number of protein queries that did not hit any hmm model ##
set num_not_found_hmm_cluster [expr $num_input_queries - $num_found_hmm_cluster]
puts $log "\tNumber not found: $num_not_found_hmm_cluster"

## Let's check the identities of the queries that did not hit anything in the bactNOG DB ##
set queries_not_hit {}
foreach query_id $query_ids {
	set presence [lsearch -glob $output_table $query_id*]
	if {$presence == -1} {
		set query_nf [lsearch -inline -glob $input_queries >$query_id*]
		lappend queries_not_hit $query_nf
	}
}

## Get the COG functional IDs from the EggNOG sqlite3 database ##
set cog_id_full_table [eggnog_cluster_to_COG_ID $output_table redundant_COG_list]
puts $log "\nNumber of hits for which COG IDs were identified: [llength $cog_id_full_table]\n\tResulting in: [llength $redundant_COG_list] redundant COG IDs"

## Write out the full table ##
set out [open $parse_process_out_name\_prot2cog.tsv w]
puts $out [join $cog_id_full_table \n]
close $out

## Write out just a list of the "narrow" COG IDs ##
set out [open $parse_process_out_name\_redundant_COGs.txt w]
puts $out [join $redundant_COG_list \n]
close $out

## For each redundant narrow entry, get the broad category ##
set broad_COG_l {}
foreach narrow_COG $redundant_COG_list {
	set broad_COG [dict get $broad_category_dict $narrow_COG]
	lappend broad_COG_l $broad_COG
}
puts [llength $redundant_COG_list]
puts [llength $broad_COG_l]
puts $log "These were translated into [llength $broad_COG_l] broad COG categories"
close $log

## Write out the list of broad categories, for plotting ##
set out [open $parse_process_out_name\_broad_COGs.txt w]
puts $out [join $broad_COG_l \n]
close $out

## Write out a list of fastas that did not get a hit against the bactNOG DB ##
set out [open $parse_process_out_name\_no_hits.faa w]
puts $out [join $queries_not_hit {}]
close $out




