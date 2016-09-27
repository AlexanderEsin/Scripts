proc parse_nog_search_output {out_file_path} {
	## A quick parser of the eggNOG output file - this parser takes only the tophit for each protein ID ##
	## Returns prot_to_nog_group_tbl: Protein_query | EggNOG_cluster_name

	## Global variables ##
	global prot_to_nog_group_tbl

	## Check that the hmmscan output file can be found ##
	if {[file exists $out_file_path] == 0} {
		puts stderr "Cannot find output file $out_file_path"
		exit 2
	}
	set output_data [split [string trim [openfile $out_file_path]] \n]

	set prot_to_nog_group_tbl {}
	set done_ids_l {}
	foreach line $output_data {
		if {[string index $line 0] eq "#"} {
			continue
		} else {
			## The output file is space-delimited, which makes it unwieldy; convert this to a tab delimited format ##
			regsub -all {\s+} $line "\t" line
			set columns [split $line \t]
			set protein_query [lindex $columns 2]

			## See if the ID has been done, this ensures that only the top hit of each result is parsed ##
			if {[lsearch $done_ids_l $protein_query] == -1} {
				## Extract the cluster name from the descriptive formatting ##
				set compound_cluster_name [lindex $columns 0]
				set processing_cluster_name [string range $compound_cluster_name [string first \. $compound_cluster_name]+1 end]
				set clean_cluster_name [string range $processing_cluster_name 0 [string first \. $processing_cluster_name]-1]

				lappend prot_to_nog_group_tbl "$protein_query\t$clean_cluster_name"
				lappend done_ids_l $protein_query
			} else {
				continue
			}
		}
	}
	return $prot_to_nog_group_tbl
}

proc eggnog_cluster_to_COG_ID {protein_2_cluster_tbl redundant_COG_list} {
	## This takes a protein -> eggNOG cluster table as produced by the parse_nog_search_output procedure, and exports two values:
	## 		1) cluster_COGcat_tbl contains: Protein_query | EggNOG_cluster_name | COG functional ID
	## 		2) redundant_COG_list: a list containing a redundant set of all the COG functional IDs in the cluster_COGcat_tbl (for plotting)
	
	## Global variables ##
	upvar $redundant_COG_list cumulative_COG_l
	global cluster_COGcat_tbl

	## Load sqlite3 package ##
	package require sqlite3

	## Load in the eggNOG database ##
	if {[file exists /users/aesin/desktop/EggNOG/eggNOG_annot_db] != 1} {
		puts stderr "Cannot find the eggNOG_annot_db Database at /users/aesin/desktop/EggNOG/"
		exit 2
	}
	sqlite3 db1 /users/aesin/desktop/EggNOG/eggNOG_annot_db

	set cluster_COGcat_tbl {}
	set cumulative_COG_l {}
	foreach hit $protein_2_cluster_tbl {
		set columns [split $hit \t]

		set protein_query [lindex $columns 0]
		set eggnogg_cluster [lindex $columns 1]

		set cog_cat [db1 eval {SELECT COGcategory from t1 where group_id = $eggnogg_cluster}]
		set function [db1 eval {SELECT function from t1 where group_id = $eggnogg_cluster}]
		
		## Check that a value is returned from the database ##
		if {[string length $cog_cat] != 0} {
			lappend cluster_COGcat_tbl "$protein_query\t$eggnogg_cluster\t$cog_cat\t$function"
		} else {
			puts stderr "EggNOG cluster: $eggnogg_cluster did NOT return a functional category. Consult the annotation table."
			exit 2
		}

		## If a joined category is returned, add both categories to the final total ##
		if {[string length $cog_cat] == 1} {
			lappend cumulative_COG_l $cog_cat
		} else {
			set cog_cats [split $cog_cat {}]
			foreach cat $cog_cats {
				lappend cumulative_COG_l $cat
			}
		}
	}
	db1 close
	return $cluster_COGcat_tbl
}

## Dictionary to "translate" the letter-based "narrow" COG IDs into broad classes ##
set broad_category_dict [dict create	\
										C Metabolism \
										G Metabolism \
										E Metabolism \
										F Metabolism \
										H Metabolism \
										I Metabolism \
										P Metabolism \
										Q Metabolism \
										D "Cellular processes and Signalling" \
										Y "Cellular processes and Signalling" \
										V "Cellular processes and Signalling" \
										T "Cellular processes and Signalling" \
										M "Cellular processes and Signalling" \
										N "Cellular processes and Signalling" \
										Z "Cellular processes and Signalling" \
										W "Cellular processes and Signalling" \
										U "Cellular processes and Signalling" \
										O "Cellular processes and Signalling" \
										J "Information Storage and Processing" \
										A "Information Storage and Processing" \
										K "Information Storage and Processing" \
										L "Information Storage and Processing" \
										B "Information Storage and Processing" \
										R "Poorly Characterized" \
										S "Poorly Characterized" \
]