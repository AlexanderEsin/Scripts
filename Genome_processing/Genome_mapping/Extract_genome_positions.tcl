#!/usr/local/bin/tclsh
source ~/Dropbox/Scripts/General_utils.tcl
source ~/Dropbox/Scripts/Procs/gc_content.tcl

package require sqlite3

## Set paths ##
set direct /users/aesin/desktop
set working_directory $direct/Geo_analysis/Geo_ortholog_nucl

set id_end_org_transl_file $direct/Clean_Proteomes/All_fastas_geo_mark_key_file.tsv

## Set folder names ##
set features_folder $direct/Clean_Genomes/Prok_feature_tables
set gbff_folder $direct/Clean_Genomes/Prok_gbff_files
set gbff_derived_prot_folder $direct/Clean_Genomes/Gbff_derived_proteomes
set gbff_derived_prot_db $direct/Clean_Genomes/Gbff_derived_proteomes_db
set fna_folder $direct/Clean_Genomes/Prok_genomic_fna_reform

set geobac_nucl_folder $working_directory/Geobac_nucl_seqs
set temp_dir $working_directory/Temp

## Open necessary databases ##
sqlite3 all_prot_geo_db $direct/Clean_Proteomes/all_prot_geo_ncbi_db
sqlite3 reannotation_report_db $direct/Clean_Proteomes/Reannotation_db
sqlite3 feature_table_db $direct/Clean_Genomes/Feature_table_db

## NEW DATABASE NAME ##
set db_new_name Geo_v_all_orth_location_NEW2_db

## Make any necessary directories ##
file mkdir $temp_dir

###########################################################################

## Create new sqlite database to hold the genomic location data ##
## Columns: Group | Prot_ID | Prot_id_no_end | ALT_ID | Accession | Start | End | Strand | Locus_tag | Binomial | Class | File_Name
cd $working_directory
if {[file exists $db_new_name] == 1} {
	set db_files [glob $db_new_name*]
	foreach file $db_files {file delete $file}
}
sqlite3 db_new $db_new_name
db_new eval {PRAGMA main.page_size = 4096}
db_new eval {PRAGMA main.cache_size=10000}
db_new eval {PRAGMA main.locking_mode=EXCLUSIVE}
db_new eval {PRAGMA main.synchronous=NORMAL}
db_new eval {PRAGMA main.journal_mode=WAL}
db_new eval {create table t1(group_number text, prot_id text, prot_id_no_end text, alt_id text, accession text, gene_start text, gene_end text, strand text, locus_tag text, binomial text, class text, file_name text, ncbi_id text)}

## Open the log file to which we will write all the results of addition ##
set log_results [open $db_new_name\_log.txt w]

## Go to test dir and pick up test file ##
cd $working_directory/Groups
set protein_families [lsort -dictionary [glob *faa]]

### COUNTERS ###
## Global ##
set total_proteins_processed 0
set total_eukaryotic_processed 0
set total_geobac_processed 0
set total_added_to_db 0
set total_geobac_added_to_db 0

set break_clause 0

set master_missing_geobac_l {}


foreach family $protein_families {

	set group_number [string range $family 0 end-4]
	puts stdout "Working on $group_number"; puts $log_results "\n######################################\nWorking on $group_number"

	
	## Open test file and get protein sequences ##
	set prot_seqs [split_genes [openfile $family]]
	set prot_seq_num [llength $prot_seqs]
	progress_init $prot_seq_num

	set multimappers {}
	set geobac_not_found_l {}


	### COUNTERS ###
	## Group specific ##
	set total_counter 0
	set added_to_db 0
	set geobac_added_to_db 0
	set no_gbff_file 0

	set euk_counter 0
	set geobac_counter 0
	set reannotated_counter 0

	set unchanged_counter 0
	set found_in_feat_counter 0
	set feat_not_present_counter 0
	set blast_miss_counter 0
	set no_blast_hit_counter 0
	set blast_hit_counter 0

	set blast_then_feat_counter 0
	set blast_locus_feat_counter 0
	set no_feat_gbk_counter 0

	set multimap_counter 0
	set no_feature_table_counter 0
	set no_feat_prot_found 0
	set no_feat_found_after_blast 0
	set no_feat_not_found_at_all 0


	foreach prot $prot_seqs {
		incr total_counter
		progress_tick $total_counter

		set prot [string trim $prot]

		set ncbi_id ""
		set accession ""
		set gene_start ""
		set gene_end ""
		set strand ""
		set locus_tag ""
		set alt_id ""
		set new_id_from_blast ""

		######################################
		## Prot ID (e.g. WP_010512312.1.123) ##
		set prot_id [string range $prot 1 [string first " " $prot]-1]

		set entry [split [all_prot_geo_db eval {SELECT class, binomial, file_name, ncbi_id FROM t1 WHERE id = $prot_id}]]
		if {[llength $entry] != 4} {
			puts stderr "Cannot find the protein_id entry in the all_prot_geo_db database. Exiting ..."
			puts $entry
			exit
		}

		## Organism name and classification (Euk, Bact, Arch or Geobac). Trim brackets from class entry. ##
		set class [string range [lindex $entry 0] 1 end-1]
		set org_name [lindex $entry 1]

		if {$class eq "Euk"} {
			#puts "Not working on eukaryotic proteins --> skipping"
			incr euk_counter
			continue
		} elseif {$class eq "Geobac"} {
			incr geobac_counter
		}

		######################################

		## ID_no_end (e.g. WP_010512312.1) ##
		set prot_id_no_end [string range $prot_id 0 [string last \. $prot_id]-1]

		## Get the file name associated with organism for id_org_transl_tbl (e.g. Acaryochloris_marina_MBIC11017.faa), without suffix (Acaryochloris_marina_MBIC11017), and use that to search the All_fastas_ncbi_id_binomial_translation table for the NCBI ID (e.g. GCF_000018105.1_ASM1810v1) ##
		set file_name [lindex $entry 2]
		set ncbi_id [lindex $entry 3]

		## If ncbi ID was not found, see if the organism name was changed ##
		if {$ncbi_id eq ""} {
			## Previously, when I used the precursor of this script to estimate GC content in genes versus their encompassing windows, I came across a problem where I had downloaded genomes using the original NCBI IDs (GCFXXX_ASMXXXX format), but their binomial name had been updated. ##
			## Thus I had a mismatch between the file name used in the original proteome and the NCBI ID. I corrected this by downloading the Taxon Name Updates since version 60 (see Taxa_name_update_tbl.tcl), and creating an "old -> new" taxon name update table. If then, I could not find the NCBI ID based on the proteome file name, I would consult this table for a change in nomenclature and use the new name to isolate the NCBI ID. ##
			## In addition, on the first run through, I created a final_name_ncbi_ID_translation_list which contained a final mapping between the proteome file name and the ncbi IDs so that I could now directly find the NCBI IDs without looking at the translation table. This new mapping was written to .../Clean_Proteomes/All_fastas_ncbi_id_binomial_translation.tsv and in turn was latter written into the all_prot_geo_ncbi_db SQLITE database ##

			## As a result, we should now never expected to not find the correct NCBI ID when querying using the protein ID ##
			## Consult the OLD_calc_gc_content_fam.tcl for examples ##

			puts stderr "NCBI ID ERROR"
			exit
		}

		set reannotated_id ""
		## Because the old IDs were always redundant, there should always be a maximum of one entry in the reannotation_db for each old id mapping to a new id. The reverse is not necessarily true. ##
		set reannotated_entry [split [reannotation_report_db eval {SELECT new_id, accession, new_start, new_end, new_strand, new_locus_tag FROM t1 WHERE old_id = $prot_id_no_end}]]
		set reannotated_id [lindex $reannotated_entry 0]

		if {$reannotated_id eq "{}"} {
			## If there is a hit, but the reannotated_id column is missing it means that in the reannotation the gene is either predicted as a pseudogene, or the entry was removed entirely. We will ignore these as the protein prouct can no longer be mapped to the genome ##
			lappend psuedogenized_list "$org_name\t$prot_id"
			incr unchanged_counter
			continue
			## SHOULD WE HAVE THIS FILTER DOWN ?? ##
			## Probably not ##
		} elseif {$reannotated_id ne ""} {
			set accession [lindex $reannotated_entry 1]
			set gene_start [lindex $reannotated_entry 2]
			set gene_end [lindex $reannotated_entry 3]
			set strand [lindex $reannotated_entry 4]
			set locus_tag [lindex $reannotated_entry 5]
			set alt_id $reannotated_id

			db_new eval {insert into t1 values($group_number,$prot_id,$prot_id_no_end,$alt_id,$accession,$gene_start,$gene_end,$strand,$locus_tag,$org_name,$class,$file_name,$ncbi_id)}
			if {$class eq "Geobac"} {
				incr geobac_added_to_db
			}
			incr reannotated_counter
			incr added_to_db
			continue
		} else {
			incr unchanged_counter
		}

		######################################
		
		## Set gbff_file and and gbff_derived_proteome_file handles ##
		set gbff_file "$gbff_folder/$ncbi_id\_genomic.gbff"
		set gbff_derived_proteome_file "$gbff_derived_prot_folder/$ncbi_id\_genomic.faa"
		## Check that the gbff file exists for this entry ##
		if {[file exists $gbff_file] == 0} {
			puts $log_results "Gbff file for $org_name --> $ncbi_id DOES NOT EXIST"
			incr no_gbff_file
			continue
		}

		######################################

		## If there is a Feature_table file, open it and look for the protein id there ##
		if {[file exists $features_folder/$ncbi_id\_feature_table.txt] == 1} {
			#puts "Searching feature table: $ncbi_id\_feature_table.txt"

			set feature_tbl_entry [split [feature_table_db eval {SELECT accession, gene_start, gene_end, strand, locus_tag FROM t1 WHERE prot_id = $prot_id_no_end AND ncbi_id = $ncbi_id}]]
			set accession [lindex $feature_tbl_entry 0]

			## Assuming the first part of the DB entry corresponds to the accession, and we pick up only one line (4 columns), we can be sure that this is the correct entry ##
			if {$accession ne "" && [llength $feature_tbl_entry] == 5} {
				set gene_start [lindex $feature_tbl_entry 1]
				set gene_end [lindex $feature_tbl_entry 2]
				set strand [lindex $feature_tbl_entry 3]
				set locus_tag [lindex $feature_tbl_entry 4]

				#puts "FEATURE_TABLE: $prot_id\t$accession\t$gene_start\t$gene_end\t$strand"
				incr found_in_feat_counter
			## If we get more than 4 columns, we must have picked up 2 or more records - this protein ID maps to several parts of the genome in this assembly ##
			} elseif {[llength $feature_tbl_entry] > 5} {
				lappend multimappers "$prot_id\t$ncbi_id"
				#puts "Incorrect length\n$prot_id_no_end\n$feature_tbl_entry"
				if {$class eq "Geobac"} {lappend geobac_not_found_l "Multimapper\t$prot_id\t$org_name\t$ncbi_id"}
				incr multimap_counter
				continue
			## Otherwise, although there is a feature table, the protein ID does not appear there (or in the reannotation DB) - this is probably a WP_XX -> WP_XX renaming ##
			} else {
				incr feat_not_present_counter

				#puts "Not found in feature table. $prot_id == $ncbi_id"
				quick_blast $prot $ncbi_id $gbff_derived_prot_db $temp_dir
				quick_process_blast $temp_dir

				if {[string length $new_id_from_blast] != 0} {
					#puts "\nBlast translated: $prot_id_no_end -> $new_id_from_blast"

					if {$new_id_from_blast eq "missing_protein_id_qualifer"} {
						#puts "Genbank file indicates a missing ID. $ncbi_id"
						incr blast_miss_counter
						if {$class eq "Geobac"} {lappend geobac_not_found_l "Missing protein ID\t$prot_id\t$org_name\t$ncbi_id"}			
						continue
					} else {
						incr blast_hit_counter
					}
					
					set feature_tbl_entry [split [feature_table_db eval {SELECT accession, gene_start, gene_end, strand, locus_tag FROM t1 WHERE prot_id = $new_id_from_blast AND ncbi_id = $ncbi_id}]]
					set accession [lindex $feature_tbl_entry 0]

					if {$accession ne "" && [llength $feature_tbl_entry] == 5} {
						#puts "\t...and it was then found in the feature table"
						set gene_start [lindex $feature_tbl_entry 1]
						set gene_end [lindex $feature_tbl_entry 2]
						set strand [lindex $feature_tbl_entry 3]
						set locus_tag [lindex $feature_tbl_entry 4]

						#puts "FEATURE_TABLE: $prot_id\t$accession\t$gene_start\t$gene_end\t$strand"
						incr blast_then_feat_counter
					## If we get more than 4 columns, we must have picked up 2 or more records - this protein ID maps to several parts of the genome in this assembly ##
					} elseif {[llength $feature_tbl_entry] > 5} {
						lappend multimappers "$prot_id\t$ncbi_id"
						#puts "Incorrect length\n$prot_id_no_end\n$feature_tbl_entry"
						if {$class eq "Geobac"} {lappend geobac_not_found_l "Multimapper\t$prot_id\t$org_name\t$ncbi_id"}		
						incr multimap_counter
						continue
					} else {
						#puts "\t...and it the protein ID was not found in feature table. Attempting to use Locus Tag..."
						## Let's find the locus tag of the blast hit in the gbff_derived_proteome ##
						set gbff_der_prot_seqs [split_genes [openfile $gbff_derived_proteome_file]]
						set gbff_der_prot_hits [lsearch -all -inline $gbff_der_prot_seqs *$new_id_from_blast*]
						
						if {[llength $gbff_der_prot_hits] == 1} {
							set gbff_der_prot_hit [lindex $gbff_der_prot_hits 0]
							set locus_tag [get_locus_tag_from_protein_header $gbff_der_prot_hit]

						} elseif {[llength $gbff_der_prot_hits] > 1} {
							#puts "\t...the protein ID ($new_id_from_blast) corresponds to multiple hits in the proteome. Multimapper"
							lappend multimappers $prot_id
							#puts "Incorrect length\n$prot_id_no_end\n$feature_tbl_entry"
							if {$class eq "Geobac"} {lappend geobac_not_found_l "Multimapper\t$prot_id\t$org_name\t$ncbi_id"}			
							incr multimap_counter
							continue

						} else {
							puts stderr "This should not happen? ERROR A1"
							puts stderr "$prot_id\t$ncbi_id\t$new_id_from_blast"
							exit
						}

						#puts "\tThe locus tag is: $locus_tag. Searching feature database..."
						## Now let's use the locus tag of the blast hit to find the relevant entry in the feature table ##
						set feature_tbl_entry [split [feature_table_db eval {SELECT prot_id, accession, gene_start, gene_end, strand FROM t1 WHERE locus_tag = $locus_tag AND ncbi_id = $ncbi_id}]]
						set accession [lindex $feature_tbl_entry 1]

						if {$accession ne "" && [llength $feature_tbl_entry] == 5} {
							#puts "\t...and it was found in the feature table"
							set gene_start [lindex $feature_tbl_entry 2]
							set gene_end [lindex $feature_tbl_entry 3]
							set strand [lindex $feature_tbl_entry 4]

							#puts "FEATURE_TABLE: $prot_id\t$accession\t$gene_start\t$gene_end\t$strand"
							incr blast_locus_feat_counter

						## If we get more than 4 columns, we must have picked up 2 or more records - this protein ID maps to several parts of the genome in this assembly ##
						} elseif {[llength $feature_tbl_entry] > 5} {
							lappend multimappers "$prot_id\t$ncbi_id"

							#puts "Incorrect length\n$prot_id_no_end\n$feature_tbl_entry"
							if {$class eq "Geobac"} {lappend geobac_not_found_l "Multimapper\t$prot_id\t$org_name\t$ncbi_id"}			
							incr multimap_counter
							continue

						## If we can't find the entry in the features table - let's ignore the existence of the feature table and just use the location provided by the genbank file. However, this will be the CDS coordinates rather than gene coordinates which we pull out from the reannotated_db or the feature tables ##
						} else {
							#puts "Feature table is not useful at all, getting location and details from genbank file ..."
							features_from_protein_header $gbff_der_prot_hit
							incr no_feat_gbk_counter
						}
					}
				} else {
					#puts "Blast did not find a hit ..."
					if {$class eq "Geobac"} {lappend geobac_not_found_l "Not found via blast\t$prot_id\t$org_name\t$ncbi_id"}			
					incr no_blast_hit_counter
				}		
			}
		## If there is not feature table, and the protein has not been reannotated, let's attempt to find the protein ID directly in the gbff-derived proteome ##
		} else {
			incr no_feature_table_counter
			## Open the fasta file, extract fasta entries, and search for the protein ID ##
			set gbff_der_prot_seqs [split_genes [openfile $gbff_derived_proteome_file]]
			set gbff_der_prot_hits [lsearch -all -inline $gbff_der_prot_seqs *$prot_id_no_end*]

			## If the protein ID matches a single entry in the proteome, we can extract the location, accession, locus_tag from the header ##
			if {[llength $gbff_der_prot_hits] == 1} {
				set gbff_der_prot_hit [lindex $gbff_der_prot_hits 0]

				set locus_tag [get_locus_tag_from_protein_header $gbff_der_prot_hit]
				features_from_protein_header $gbff_der_prot_hit
				incr no_feat_prot_found

			## If we get more than one hit, two genes express the same protein product and it's not going to be useful for us ##
			} elseif {[llength $gbff_der_prot_hits] > 1} {
				#puts "\t...the protein ID ($prot_id_no_end) corresponds to multiple hits in the proteome. Multimapper"
				lappend multimappers $prot_id
				#puts "Incorrect length\n$prot_id_no_end\n$feature_tbl_entry"
				if {$class eq "Geobac"} {lappend geobac_not_found_l "Multimapper\t$prot_id\t$org_name\t$ncbi_id"}			
				incr multimap_counter
				continue

			## If we do not find it, perhaps the protein ID changed, so let's try blasting ##
			} else {

				quick_blast $prot $ncbi_id $gbff_derived_prot_db $temp_dir
				quick_process_blast $temp_dir

				if {[string length $new_id_from_blast] != 0} {
					#puts "\nBlast translated: $prot_id_no_end -> $new_id_from_blast"

					if {$new_id_from_blast eq "missing_protein_id_qualifer"} {
						#puts "Genbank file indicates a missing ID. $ncbi_id"
						# incr blast_miss_counter
						if {$class eq "Geobac"} {lappend geobac_not_found_l "Missing protein ID\t$prot_id\t$org_name\t$ncbi_id"}			
						continue
					} else {
						# incr blast_hit_counter
					}

					set gbff_der_prot_hits [lsearch -all -inline $gbff_der_prot_seqs *$new_id_from_blast*]
					
					if {[llength $gbff_der_prot_hits] == 1} {
						set gbff_der_prot_hit [lindex $gbff_der_prot_hits 0]

						set locus_tag [get_locus_tag_from_protein_header $gbff_der_prot_hit]
						set alt_id $new_id_from_blast
						features_from_protein_header $gbff_der_prot_hit
						incr no_feat_found_after_blast

					} elseif {[llength $gbff_der_prot_hits] > 1} {
						#puts "\t...the protein ID ($new_id_from_blast) corresponds to multiple hits in the proteome. Multimapper"
						lappend multimappers $prot_id
						#puts "Incorrect length\n$prot_id_no_end\n$feature_tbl_entry"
						if {$class eq "Geobac"} {lappend geobac_not_found_l "Multimapper\t$prot_id\t$org_name\t$ncbi_id"}	
						incr multimap_counter
						continue
					} else {
						puts stderr "This should not happen? ERROR A3"
						puts stderr "$prot_id\t$ncbi_id\t$new_id_from_blast"
						exit
					}
				} else {
					if {$class eq "Geobac"} {lappend geobac_not_found_l "Not found via blast\t$prot_id\t$org_name\t$ncbi_id"}	
					incr no_feat_not_found_at_all
				}
			}
		}

		## Remove any double dashes ##
		if {[regexp "__" $org_name] == 1 || [regexp "__" $file_name] == 1} {
			puts "Double dashes present: $file_name\t$org_name"
			regsub -all "__" $org_name "_" $org_name
			regsub -all "__" $file_name "_" $file_name
			puts "Double dash replacement: $file_name\t$org_name"
		}

		if {$strand ne "+" && $strand ne "-" && $strand ne "" && $strand ne "?"} {
			puts "STRAND ERROR: $strand"
			set break_clause 1
			break
		}

		if {$accession eq "" || $gene_start eq "" || $gene_end eq "" || $strand eq ""|| $locus_tag eq ""} {
			incr not_addeed_to_db
		} else {
			db_new eval {insert into t1 values($group_number,$prot_id,$prot_id_no_end,$alt_id,$accession,$gene_start,$gene_end,$strand,$locus_tag,$org_name,$class,$file_name,$ncbi_id)}
			incr added_to_db
			if {$class eq "Geobac"} {
				incr geobac_added_to_db
			}
		}

	#puts "Done: $total_counter / $prot_seq_num"
	}
	if {$break_clause == 1} {break}

	set master_missing_geobac_l [concat $master_missing_geobac_l $geobac_not_found_l]
	puts stdout "Cumulative num. Geobacillus proteins no coodinates: [llength $master_missing_geobac_l]"

	puts $log_results "Total proteins processed: $total_counter"
	puts $log_results "Total geobac processed: $geobac_counter"
	puts $log_results "Do not have gbff files: $no_gbff_file"
	puts $log_results "Total added to DB: $added_to_db"
	puts $log_results "Total geobac added to DB: $geobac_added_to_db"
	puts $log_results "Eukaryotic proteins: $euk_counter"

	puts $log_results "Reannotated: $reannotated_counter"

	puts $log_results "Not reannotated: $unchanged_counter"
	puts $log_results "\tFound successfully in feature_db: $found_in_feat_counter"

	puts $log_results "\tNot found in feature_db: $feat_not_present_counter"
	puts $log_results "\t\tHad a blast-hit to a missing ID: $blast_miss_counter"
	puts $log_results "\t\tHad no blast-hit at all: $no_blast_hit_counter"
	puts $log_results "\t\tHad a blast-hit to a protein-ID: $blast_hit_counter"

	puts $log_results "\t\t\tWas then found in feature_table: $blast_then_feat_counter"
	puts $log_results "\t\t\tHad to use locus tag from blast hit: $blast_locus_feat_counter"
	puts $log_results "\t\t\tNot found even with locus, used genbank file: $no_feat_gbk_counter"

	puts $log_results "No feature table for these entries: $no_feature_table_counter"
	puts $log_results "\tProtein ID found in the gbff-derived proteome: $no_feat_prot_found"
	puts $log_results "\tProtein ID found in the gbff-derived proteome AFTER blast: $no_feat_found_after_blast"
	puts $log_results "\tProtein not found at all: $no_feat_not_found_at_all"

	puts $log_results "\tMap to multiple genomic loci: $multimap_counter"
	puts $log_results "Geobacilli missing:\n[join $geobac_not_found_l \n]"
	chan flush $log_results

	set total_proteins_processed [expr $total_proteins_processed + $total_counter]
	set total_geobac_processed [expr $total_geobac_processed + $geobac_counter]
	set total_eukaryotic_processed [expr $total_eukaryotic_processed + $euk_counter]

	set total_added_to_db [expr $total_added_to_db + $added_to_db]
	set total_geobac_added_to_db [expr $total_geobac_added_to_db + $geobac_added_to_db]

}
close $log_results

puts stdout "\nTotal proteins processed: $total_proteins_processed"
puts stdout "Total Geobac proteins processed: $total_geobac_processed"
puts stdout "Total eukaryotic proteines processed: $total_eukaryotic_processed"
puts stdout "\nTotal added to database: $total_added_to_db"
puts stdout "Total Geobac added to database: $total_geobac_added_to_db"

db_new eval {CREATE INDEX prot_id_index ON t1 (prot_id)}
db_new eval {CREATE INDEX locus_tag_index ON t1 (group_number)}
db_new eval {CREATE INDEX group_index ON t1 (group_number)}

#puts "\n[join $psuedogenized_list \n]"
#puts "\n[join $multimappers \n]"

all_prot_geo_db close
reannotation_report_db close
feature_table_db close
db_new close

# #################################################################################################################################################################################################################################
# #################################################################################################################################################################################################################################
# #################################################################################################################################################################################################################################

