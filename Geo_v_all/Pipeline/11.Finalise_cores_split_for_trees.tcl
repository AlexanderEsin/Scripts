#!/usr/local/bin/tclsh
###########################################################################
### Changelog ###
# 15 January 2016 - Complete revamp. Increase speed and simplicity. Fixed variable name errors in amino acid replacing #
# 29 January 2016 - Added output a list of dud sequences (no overlap with aligned core Geobacillus sequence) for inspection #


###########################################################################
### Procedures and packages ###
source ~/Dropbox/Scripts/General_utils.tcl

###########################################################################
### Set variables ###
set ival 2.0
# set direct /users/aesin/desktop/Geo_analysis/Geo_v_all/$ival
set direct /users/aesin/desktop/Geo_analysis/Geo_v_all/$ival/Missing_geobac_reconstruction

set aa_swap_switch 1
set align_pat_switch 0
set MSA_split_switch 0
set all_in_one_switch 1; # This doesn't require either align_pat_switch or MSA_split_switch

###########################################################################
### For those files that had a core alignmet, use the 'core sequence' (the sequence on which the core alignment was made) to trim down the MSA to the aligned core sequence + 25%. This will trim long tails from poorly aligned sequences, and so reduce the resolution of distantly related sequences, but significantly improve tree-building time ###
set input_path "$direct/Group_fastas_MSA"
set all_in_one_dir $direct/MSA_split

##############################
### Make output directories ###
file mkdir $input_path/Processed_UNTRIM
file mkdir $input_path/Processed_win_id

## Counters ##
set done_counter 1
set dud_counter 0
set geo_dud_counter 0

## Patterns for regexp searches ##
set pata {(k)+?.+?}

## Lists ##
set dud_list {}

## Get input files ##
cd $input_path
set untrimmed [glob -nocomplain UNTRIM*]

if {[llength $untrimmed] > 0} {

	foreach file $untrimmed {

		set output_name [string range $file 7 end]
		set dir [string range $output_name 0 end-4]
		
		## Pull out the relevant data from the win_id file and find the kX (k_value) number of the core sequence in alignment ##
		set win_id_data [split [string trim [openfile win_id_$dir]] \n]

		set win_id [string trim [lindex [split [lindex $win_id_data 0] \t] 1]]
		set win_length [string trim [lindex [split [lindex $win_id_data 1] \t] 1]]

		## Find the correct k-value for the winning ID ##
		set key_data [openfile KEY_$dir\.txt]
		set key_list [split $key_data \n]
		regexp -line $pata$win_id $key_data key_line_match
		set win_k_value [string trim [lindex $key_line_match 0]]

		##############################
		## Get the CDS of the core sequence and find the indices of the first and last residues ##
		set alignment_data [openfile $file]
		set aligned_seqs [split_genes_fast $alignment_data]

		set core_gene [lsearch -inline -glob $aligned_seqs >$win_k_value\n*]
		regsub -all "\n" [string trim [string range $core_gene [string first \n $core_gene] end]] {} CDS
		set max_length [string length $CDS]

		## Index of first amino acid ##
		regexp {[A-Z]} $CDS firstletter
		set CDS_start [string first $firstletter $CDS]

		## Index of last amino acid##
		regexp {[A-Z]} [string reverse $CDS] lastletter
		set CDS_end [string last $lastletter $CDS]
		
		## Alignment sequence between the first and last amino acids ##
		set core_CDS [string range $CDS $CDS_start $CDS_end]

		##############################
		## For each aligned sequence in the alignment, reduce the size down to core sequence +/- 25% length. This will often result in a longer "core" than the original as the core sequence in the final alignment can contain gaps. By using win_length here to calculate the +/-25% - we are not capturing the +/- 25% of the core alignment, but less because of the gaps ##

		## Set the new start and end cutoffs ##
		set core_start [::tcl::mathfunc::round [expr $CDS_start-($win_length*0.25)]]
		if {$core_start < 0} {
			set core_start 0
		}
		set core_end [::tcl::mathfunc::round [expr $CDS_end+($win_length*0.25)]]
		if {$core_end > $max_length} {
			set core_end $max_length
		}

		set new_alignment {}

		foreach seq $aligned_seqs {
			set firstn [string first \n $seq]
			set header [string trim [string range $seq 0 $firstn]]

			regsub -all "\n" [string trim [string range $seq $firstn end]] {} CDS
			set new_CDS [string range $CDS $core_start $core_end]

			## If the new aligned sequences (the core region) contains only gaps - i.e. the sequence is really diverged from the geobac of interest - do not include the sequence in the alignment ##
			if {[regexp -all -- {-} $new_CDS] != [string length $new_CDS]} {
				set line_broken_CDS [linebreak_align $new_CDS]
				lappend new_alignment "$header\n[string trim $line_broken_CDS]"

			## If there are only gaps in the aligned sequence, we add it to the dud list - i.e. sequences that share no overlap with the aligned core Geobacillus sequence. If these are also Geobacillus sequences, we extract them into a separate file for inspection ##
			} else {
				set k_val [string range $header 1 end]
				set dud_entry [split [lsearch -inline -glob $key_list $k_val\t*] \t]
				
				if {[string length $dud_entry] == 0} {puts "Can't find entry for $kval in $file"; continue}

				set dud_class [lindex $dud_entry 2]

				if {$dud_class eq "Geobac"} {
					set dud_id [lindex $dud_entry 3]
					lappend dud_list "$dir\t$dud_id"
					incr geo_dud_counter
				}

				incr dud_counter
			}
		}

		##############################
		## Removing them from the key was annoying ##

		### If any proteins were poorly aligned and were removed in the above step (all gaps over the region of the trimmed alignment), they need to be removed from the respective key. 71215 - removing the entries from the key means that we cannot track back "lost" sequence. Might as well leave them in. ###

		# set key_entries [split [string trim $key_data] \n]
		# # List length of the trimmed alignment versus the length of the key #
		# if {[llength $new_alignment] != [llength $key_entries]} {
		# 	foreach entry $key_entries {
		# 		set k_val [lindex [split $entry \t] 0]
		# 		# If the k_val in the key is not found in the alignment, remove that entry from the alignment
		# 		if {[lsearch -glob $new_alignment >$k_val\n*] == -1} {
		# 			set key_entries [lsearch -all -inline -not -exact $key_entries $entry]
		# 		}
		# 	}
		# 	set update_key_data [join $key_entries \n]
		# 	set out [open KEY_$dir\.txt w]
		# 	puts $out $update_key_data
		# 	close $out
		# }
		
		##############################
		## Output the new alignment ##

		set out [open fasta_$output_name w]
		puts $out [join $new_alignment \n]
		close $out

		set out [open $direct/Dud_list.tsv w]
		puts $out [join $dud_list \n]
		close $out

		## Convert it to phylip format for RAxML and delete the fasta format file (space) ##
		catch {exec seqret fasta_$output_name $output_name -osformat phylip}
		file delete fasta_$output_name
		
		## Stash away the "untrimmed" alignment and win_id files ##
		file rename -force $file Processed_UNTRIM/$file
		file rename -force win_id_$dir Processed_win_id/win_id_$dir

		puts "Trimming UNTRIMMED files: $done_counter / [llength $untrimmed]"
		incr done_counter
	}
}

puts "Done cutting down alignments.\nThere were $dud_counter sequences that had no overlap with the aligned 'core' Geobacillus sequences.\nOf these 'duds', $geo_dud_counter were Geobacillus sequences"


###########################################################################
### Convert all selenocysteine (U) into cysteine (C) and other non-standard AAs into their most similar standard counterparts. Also sort away those alignments containing fewer than four sequences - e.g. if an alignment had 4 sequences, and then one was removed in the alignment trimming above ###

if {$aa_swap_switch == 1} {
	cd $direct/Group_fastas_MSA
	file mkdir Insufficient_taxa

	set i 1
	set too_small_counter 0
	set fas_files [glob *.fas]

	foreach fas_file $fas_files {
		set sequence_data [openfile $fas_file]
		set x 0

		## Search for all the non-standard residues and convert them to the closest amino acid ##
		if {[regexp -all "U" $sequence_data] > 0} {
			set x [regsub -all "U" $sequence_data "C" sequence_data]
			puts "$fas_file === No. selenocysteine subs: $x"
		}
		if {[regexp -all "J" $sequence_data] > 0} {
			set x [regsub -all "J" $sequence_data "I" sequence_data]
			puts "$fas_file === No. J residues subs: $x"
		}
		if {[regexp -all "B" $sequence_data] > 0} {
			set x [regsub -all "B" $sequence_data "D" sequence_data]
			puts "$fas_file === No. B residues subs: $x"
		}
		if {[regexp -all "Z" $sequence_data] > 0} {
			set x [regsub -all "Z" $sequence_data "E" sequence_data]
			puts "$fas_file === No. Z residues subs: $x"
		} 
		if {[regexp -all "O" $sequence_data] > 0} {
			set x [regsub -all "O" $sequence_data "K" sequence_data]
			puts "$fas_file === No. O residues subs: $x"
		}

		## If any changes were made, re-write out ##
		if {$x > 0} {
			set out [open $fas_file w]
			puts $out $sequence_data
			close $out
		}

		## Filter away any alignments with fewer than 4 sequences. NB searching in phylip format, hence regexp ##
		if {[regexp -all -line {k+?.+?\n+?} $sequence_data] < 4} {
			file rename $fas_file Insufficient_taxa/$fas_file
			incr too_small_counter
		}
		puts "$i / [llength $fas_files]"
		incr i
	}
	puts "$too_small_counter groups had fewer than 4 sequences (insufficient for tree-building)"
}


###########################################################################
### Process each RAxmL-eligible alignment and collect the number of 'distinct alignment patterns' from the RAxML info file. Put all these values into the Distinct_align_patterns_2.5.tsv file. ###

if {$align_pat_switch == 1} {
	cd $direct/Group_fastas_MSA

	set aligns [lsort -dictionary [glob *fas]]
	set pattern_map 0

	if {[file exists $direct/Distinct_align_patterns.tsv] == 1} {
		openfile $direct/Distinct_align_patterns.tsv
		set pattern_map [split [string trim $data] \n]
	}

	if {[llength $aligns] != [llength $pattern_map]} {
		if {[file exists $direct/Distinct_align_patterns.tsv] == 1} {
			file delete $direct/Distinct_align_patterns.tsv
		}

		foreach align $aligns {
			regsub {.fas} $align {} out_number
			set pid [exec -keepnewline raxml -f a -s $align -n $out_number\.txt -m PROTCATDAYHOFF -p 1234 -x 10001 -N 1000 -T 4 &]
			set info [glob -nocomplain *info*]
			while {[llength $info] < 1} {
				after 100
				set info [glob -nocomplain *info*]
			}

			set info [lindex [glob -nocomplain *info*] 0]
			openfile $info
			set info_done 0
			while {$info_done < 1} {
				after 100
				openfile $info
				set info_done [string first "RAxML was called as follows" $data]
			}

			catch {exec kill $pid}

			set pattern_no [string trim [string range $data [string first "Alignment has" $data] [string first " distinct alignment" $data]]]
			set pattern_no [string trim [string range $pattern_no [string last " " $pattern_no] end]]

			set out [open $direct/Distinct_align_patterns.tsv a]
			puts $out "$out_number\t$pattern_no"
			close $out
			
			set reduced [lindex [glob -nocomplain *reduced] 0]
			file delete $reduced

			set raxmls [glob RAx*]
			foreach raxml $raxmls {
				file delete $raxml
			}
		}
	}
}

###########################################################################
### Package the MSA and keys into directories ready for tree-building ###

if {$MSA_split_switch == 1} {
	cd $direct; file mkdir MSA_split
	cd $direct/MSA_split
	file mkdir individual
	file mkdir 8_core/small; file mkdir 8_core/large; file mkdir 8_core/giant; file mkdir 8_core/colossal
	file mkdir 4_core/small; file mkdir 4_core/large; file mkdir 4_core/giant; file mkdir 4_core/colossal

	set 4_core {}; set 8_core {}; set individual {}

	cd $direct
	openfile $direct/Distinct_align_patterns.tsv
	set data [split [string trim $data] \n]
	foreach align $data {
		set align [split [string trim $align] \t]
		set group_no "[lindex $align 0].fas"
		set align_pat [lindex $align 1]
		if {$align_pat > 4000} {
			lappend individual $group_no
		} elseif {$align_pat > 1250} {
			lappend 8_core $group_no
		} else {
			lappend 4_core $group_no
		}
	}

	## SORT the biggest align_pat alignments into a folder to be done manually ##
	cd $direct/Group_fastas_MSA
	foreach file $individual {
		regsub {.fas} $file {} number
		file copy -force $file $direct/MSA_split/individual/
		file copy -force KEY_$number\.txt $direct/MSA_split/individual/
	}

	## SUB-SORT the rest of the files based on the number of sequences ##
	####
	set 4_core_parallel 64; #Number of parallel implementations
	set 8_core_parallel 12; #Number of parallel implementations
	####

	set core_distribs {4_core 8_core}

	foreach group $core_distribs {
		set x 1
		foreach file [expr $$group] {
			regsub {.fas} $file {} number
			openfile $direct/Group_fastas_MSA/$file
			set data [string trim $data]; set seq_no [string trim [string range $data 0 [string first " " $data]]]
			if {$seq_no > 3000} {
				set size colossal
			} elseif {$seq_no >= 1250} {
				set size giant
			} elseif {$seq_no >= 750} {
				set size large
			} else {
				set size small
			}

			# Make the necessary subfolders for parallelisation - if needed #
			if {[llength [glob -nocomplain -type d $direct/MSA_split/$group/$size/*]] == 0} {
				set i 1
				while {$i <= [expr $$group\_parallel]} {
					file mkdir $direct/MSA_split/$group/$size/$i
					incr i
				}
			}
			# Sort the files to the folder containing the least files #
			set input_folders [glob -type d $direct/MSA_split/$group/$size/*]
			set smallest_folder ""
			set smallest_folder_len "foo"
			foreach folder $input_folders {
				cd $folder
				set len_folder [llength [glob -nocomplain -type d *]]
				if {$len_folder < $smallest_folder_len} {
					set smallest_folder $folder
					set smallest_folder_len $len_folder
				}
				cd ..
			}
			file mkdir $smallest_folder/$number
			cd $direct/Group_fastas_MSA
			file copy -force $file $smallest_folder/$number/
			file copy -force KEY_$number\.txt $smallest_folder/$number/

			puts "$group == $x / [llength [expr $$group]] --> $size"
			incr x
		}
	}
}

###########################################################################
if {$all_in_one_switch == 1} {
	set counter 0
	file mkdir $all_in_one_dir

	cd $direct/Group_fastas_MSA
	set alignments [glob *fas]

	foreach alignment $alignments {
		set group_number [string range $alignment 0 end-4]
		file mkdir $all_in_one_dir/$group_number

		file copy -force $alignment $all_in_one_dir/$group_number
		file copy -force KEY_$group_number\.txt $all_in_one_dir/$group_number

	}
}


###########################################################################

# cd $direct
# catch {exec gzip -r MSA_split}
# catch {exec scp -r MSA_split ade110@ax3.hpc.ic.ac.uk:/scratch/ade110/Geo_v_all}
# puts "DONE"

###########################################################################