#!/usr/local/bin/tclsh
##
puts [exec pwd]
##
set direct /users/aesin/desktop/Consensus_trees
set org Geobac/Outgroups/With_excl_geos
set eval -10
set ival 2
set rbbh_dir Rbbh_all_$eval
###

proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}

proc split_genes {fasta} {
	global genes
	regsub -all {>} [string trim $fasta] {£>} fasta
	set fasta [split $fasta £]
	regsub -all "{}" $fasta {} genes
}

package require sqlite3
package require math::statistics

###########################################################################
file mkdir $direct/$org/$eval
### Make a filtered glist here from the mcl groups ###
cd $direct/$org/$eval
if {[file exists Group_fastas_$ival] == 1} {
	file delete -force Group_fastas_$ival
}
file mkdir Group_fastas_$ival

### Open the MCL-output groups file ###
cd $direct/$org/$rbbh_dir
openfile Master_groups_weight_$ival\.txt
set glist $data

### Open the prot database & create Group_fastas output folder ###
cd /users/aesin/desktop/Clean_proteomes
sqlite3 db1 all_prot_geo_ncbi_db
###

set g 1
### Identify a list of the expected binomials by sampling from the sqlite database using the Input/Fastas file list as a guide ###
set full_binomial_set {}
cd $direct/$org/Input/Fastas
if {[llength [glob -nocomplain *gz]] > 0} {
	cd $direct/$org/Input
	catch {exec gunzip -r Fastas}
}
set input_files [glob -nocomplain *.faa]
set org_number [llength $input_files]
if {[llength $input_files] == 0} {
	puts "NEED SOME INPUT FASTA FILES"
	puts "Current directory: [exec [pwd]]"
	error 2
}

foreach file $input_files {
	db1 eval {SELECT binomial FROM t1 WHERE file_name = $file} {
		set binomial "$binomial"
	}
	lappend full_binomial_set $binomial
}
puts $full_binomial_set
puts [llength $full_binomial_set]

### For each group (gene family) pull out the class (Bact/Arch etc) and full gene to create a (filterable) gene list (Group_list.tsv) and fastas for each group (Group_fastas) ###
set genefams [split $glist \n]
regsub -all "{}" $genefams {} genefams
foreach fam $genefams {
	set genes [split $fam \t]
	## REMEMBER TO CHANGE TO NUMBER OF TOTAL ORGANISMS ##
	if {[llength $genes] < $org_number} {
		puts "Too short"
	} else {
		set i 1
		set desig_fam ""
		foreach gene $genes {
			db1 eval {SELECT class FROM t1 WHERE id = $gene} {
				set class "$class"
			}
			append desig_fam "$gene $class\t"
			db1 eval {SELECT gene FROM t1 WHERE id = $gene} {
				set new_gene "$gene"
			}
			set out [open $direct/$org/$eval/Group_fastas_$ival/$g.faa a]
			puts $out $new_gene
			close $out
			puts "Group: $g/[llength $genefams] === $i/[llength $genes]"
			incr i
		}
		set out [open $direct/$org/$eval/Group_list_$ival.tsv a]
		puts $out "$g\t$desig_fam"
		close $out
	}
	incr g
}

###########################################################################
### Convert the binomial names for each gene into a phylip-friendly format + fix MULTISPECIES name tags ###
### All the regsubs attempt to remove random crap from the binomial name to avoid complications downstream ###

cd $direct/$org/$eval
file mkdir Group_fastas_key_ids_$ival
file mkdir Group_fastas_paralogs_$ival

cd $direct/$org/$eval/Group_fastas_$ival

set fams [glob *.faa]
set fams [lsort -dictionary $fams]
set i 1

set delete_list ""
set missing_list ""

foreach fam $fams {
	set missing_break_flag 0
	set paralog_break_flag 0
	set fam_binomials {}
	set fam_paralog_list {}
	set paralog_counter 0
	set present_taxa {}
	set missing_taxa {}

	openfile $fam
	set new_glist {}
	split_genes $data
	foreach gene $genes {
		set gene [string trim $gene]
		regsub -all {\*} $gene {} gene
		set binomial [string range [string range $gene 0 [string first \n $gene]-1] [string last \[ $gene]+1 end-1]
		### If there is a MULTISPECIES tag, then query the all_prot_geo_db with the gene_id to get the binomial of the organism ###
		if {[regexp {MULTISPECIES} $gene] == 1} {
			regexp {>+?.+?\s+?} $gene query
			set query [string trim [string range $query 1 end]]
			db1 eval {SELECT binomial FROM t1 WHERE id = $query} {
				set binomial "$binomial"
			}
			# The [binomial name] occasionally contains (bracketed terms) #
			regsub -all {\(} $binomial {} binomial
			regsub -all {\)} $binomial {} binomial
			regsub -all {\:} $binomial {} binomial
			regsub -all {\#} $binomial {} binomial
			regsub -all {\=} $binomial {} binomial
			regsub -all {\/} $binomial {} binomial
			regsub -all {'} $binomial {} binomial
			regsub -all " " $binomial {_} match
			regsub -all {\.} $match {} match
			regsub -all {\,} $match {} match
			regsub -all {\+} $match {} match
			regsub -all -- {\-} $match {_} match
			#
			set header_trim [string range $gene 0 [string last \[ $gene]-1]
			set new_header "$header_trim\[$match\]"
			set gene2 "$new_header\n[string trim [string range $gene [string first \n $gene] end]]"
			lappend new_glist $gene2
		} else {
			# The [binomial name] occasionally contains (bracketed terms) #
			regsub -all {\(} $binomial {} binomial
			regsub -all {\)} $binomial {} binomial
			regsub -all {\:} $binomial {} binomial
			regsub -all {\#} $binomial {} binomial
			regsub -all {\=} $binomial {} binomial
			regsub -all {\/} $binomial {} binomial
			regsub -all {'} $binomial {} binomial
			regsub -all " " $binomial {_} match
			regsub -all {\.} $match {} match
			regsub -all {\,} $match {} match
			regsub -all {\+} $match {} match
			regsub -all -- {\-} $match {_} match
			set header_trim [string range $gene 0 [string last \[ $gene]-1]
			set new_header "$header_trim\[$match\]"
			set gene2 "$new_header\n[string trim [string range $gene [string first \n $gene] end]]"
			lappend new_glist $gene2
		}
		## Find families where there is more than one copy of the same binomial ##
		if {[lsearch $fam_binomials $match] == -1} {
			lappend fam_binomials $match
		} else {
			set paralog_break_flag 1
			if {[lsearch $fam_paralog_list $match] == -1} {
				lappend fam_paralog_list $match
			}
			incr paralog_counter
		}
		## Every taxon in the family is added to a "present taxa" list ##
		lappend present_taxa $match
	}
	### Check all the present taxa in the family against the full set of 'expected' taxa - and identify those which are missing ###
	foreach taxon $full_binomial_set {
		if {[lsearch $present_taxa $taxon] == -1} {
			lappend missing_taxa $taxon
			set missing_break_flag 1
		}
	}
	set missing_taxa [join $missing_taxa "\t"]
	append missing_list "$fam\t$missing_taxa\n"

	###
	set total_unique_proteins [expr [llength $genes] - $paralog_counter]
	if {$total_unique_proteins == $org_number} {
		## PUT SOMETHING HERE? ##
	}
	### Output the new fasta file to a temp directory (originally for debugging) ###
	if {$missing_break_flag == 0 && $paralog_break_flag == 0} {
		set new_glist_string [string trim [join $new_glist \n]]
		set out [open $direct/$org/$eval/Group_fastas_key_ids_$ival/$i.faa w]
		puts $out $new_glist_string
		close $out
	} elseif {$missing_break_flag == 0 && $paralog_break_flag == 1} {
		set new_glist_string [string trim [join $new_glist \n]]
		set paralog_out ""
		set paralog_seqs {}
		split_genes $new_glist_string
		foreach paralog $fam_paralog_list {
			foreach gene $genes {
				if {[string match *$paralog* $gene] == 1} {
					lappend paralog_seqs $gene
				}
			}
			lappend paralog_seqs "\n"
		}
		set paralog_out [string trim [join $paralog_seqs {}]]
		set out [open $direct/$org/$eval/Group_fastas_paralogs_$ival/$i.faa w]
		puts $out $paralog_out
		close $out

		### Try to process paralogs to include more sequences. E.g. if one of the paralogs is less than half the length of the other pick the longer one. Else if there is a "partial" identifier or one is MULTISPECIES, pick the other sequence. ###

		regsub -all "\n\n" $paralog_out {&} paralog_out_easy_split
		set paralog_pairs [split $paralog_out_easy_split &]
		set pair [lindex $paralog_pairs 0]
		foreach pair $paralog_pairs {
			###################
			### Process the paralogs to try and filter out poor sequences ###
			set pair [string trim $pair]
			split_genes $pair
			if {[llength $genes] != 2} {
				continue
			}
			set processing_success_counter 0
			while {$processing_success_counter == 0} {
				### process partials ###
				set partials [lsearch -all -glob $genes *partial*]
				if {[llength $partials] == 1} {
					set bad_gene [string trim [lindex $genes $partials]]
					set bad_gene_id [string trim [string range $bad_gene 0 [string first " " $bad_gene]]]
					set bad_gene_index [lsearch -glob $new_glist $bad_gene_id*]
					set new_glist [lreplace $new_glist $bad_gene_index $bad_gene_index]
					
					set processing_success_counter 1
					continue
				}
				### process MULTISPECIES ###
				set multispecs [lsearch -all -glob $genes *MULTISPECIES*]
				if {[llength $multispecs] == 1} {
					set bad_gene [string trim [lindex $genes $multispecs]]
					set bad_gene_id [string trim [string range $bad_gene 0 [string first " " $bad_gene]]]
					set bad_gene_index [lsearch -glob $new_glist $bad_gene_id*]
					set new_glist [lreplace $new_glist $bad_gene_index $bad_gene_index]

					set processing_success_counter 1
					continue
				}
				### process very short sequnce stubs ###
				set paralog_lens {}
				foreach gene $genes {
					set id [string trim [string range $gene 0 [string first " " $gene]]]
					set CDS [string trim [string range $gene [string first \n $gene] end]]
					regsub -all "\n" $CDS {} CDS_flat
					lappend paralog_lens "[string length $CDS_flat]\t$id"
				}
				set paralog_lens [lsort $paralog_lens]
				set longest_paralog [string trim [lindex [lindex $paralog_lens 0] 0]]
				set 50_pc_cutoff [::tcl::mathfunc::round [expr $longest_paralog / 2]]
				if {[string trim [lindex [lindex $paralog_lens 1] 0]] < $50_pc_cutoff} {
					set bad_gene_id [string trim [lindex [lindex $paralog_lens 1] 1]]
					set bad_gene_index [lsearch -glob $new_glist $bad_gene_id*]
					set bad_gene [lindex $new_glist $bad_gene_index]
					set new_glist [lreplace $new_glist $bad_gene_index $bad_gene_index]

					set processing_success_counter 1
					continue
				}

				if {$processing_success_counter == 0} {
					set processing_success_counter 2
				}
			}
			
			if {$processing_success_counter == 1} {
				openfile $direct/$org/$eval/Group_fastas_paralogs_$ival/$i.faa
				set data "[string trim $data]\n\n###Gene Removed###\n$bad_gene"
				set out [open $direct/$org/$eval/Group_fastas_paralogs_$ival/$i.faa w]
				puts $out $data
				close $out
			}
		}
		###################

		if {[llength $new_glist] == [llength $full_binomial_set]} {
			set new_glist_string [string trim [join $new_glist \n]]
			set out [open $direct/$org/$eval/Group_fastas_key_ids_$ival/$i.faa w]
			puts $out $new_glist_string
			close $out
		}

		###
		set temp_delete_list [join $fam_paralog_list "\t"]
		append delete_list "$fam\t$temp_delete_list\n"
	}
	puts "Done: $i/[llength $fams]"
	incr i
}

set out [open $direct/$org/$eval/Ortholog_groups_not1to1_$ival.txt w]
puts $out [string trim $delete_list]
close $out
set out [open $direct/$org/$eval/Missing_taxa_$ival.txt w]
puts $out [string trim $missing_list]
close $out

###########################################################################

cd $direct/$org/$eval/Group_fastas_key_ids_$ival
set ortho_files [lsort -dictionary [glob *faa]]
set master_output "File Name\tProtein Name\tOver 10% 'hypothetical'?\tMean length\tStandard dev\tNo. len >/< 2*stdev\tNo. len >/< 8*stdev\tList of >8stdev genes\n"

foreach file $ortho_files {
	openfile $file
	set data [string trim $data]

	#### Get the protein family name and calculate how many of these are "hypothetical" proteins ####
	set first_line [string range $data 0 [string first \n $data]]
	set prot_name [string range $first_line [string first \) $first_line] [string first \[ $first_line]]
	set prot_name [string trim [string range $prot_name 1 end-1]]
	set hypotheticals [regexp -all {hypothetical} $data]
	set trim_data $data
	while {[string match -nocase *hypothetical* $prot_name] > 0 || [string match -nocase "*membrane protein*" $prot_name] > 0} {
		set trim_data [string range $trim_data 1 end]
		set trim_data [string trim [string range $trim_data [string first > $trim_data] end]]
		set first_line [string range $trim_data 0 [string first \n $trim_data]]
		set prot_name [string range $first_line [string first \) $first_line] [string first \[ $first_line]]
		set prot_name [string trim [string range $prot_name 1 end-1]]
	}
	set output_line "$file\t$prot_name"
	if {$hypotheticals > [expr $org_number / 10]} {
		append output_line "\tYES"
	} else {
		append output_line "\tNO"
	}

	#### Calculate the standard deviation of the protein length and see whether there are many outliers ####

	set aa_lengths {}
	set id_len_dict {}
	split_genes $data
	foreach gene $genes {
		set gene [string trim $gene]
		set gene_id [string trim [string range $gene 1 [string first " " $gene]]]
		set CDS [string trim [string range $gene [string first \n $gene] end]]
		regsub -all "\n" $CDS {} CDS
		set aa_length [string length $CDS]
		lappend aa_lengths $aa_length
		lappend id_len_dict "$gene_id\t$aa_length"
	}
	set stdev [::tcl::mathfunc::round [::math::statistics::stdev $aa_lengths]]
	set mean [::tcl::mathfunc::round [::math::statistics::mean $aa_lengths]]

	set over_2stdev {}
	set over_8stdev {}
	set over_8stdev_list ""
	foreach gene $id_len_dict {
		set len [lindex $gene 1]
		if {$len > [expr $mean+(2*$stdev)] || $len < [expr $mean-(2*$stdev)]} {
			lappend over_2stdev [lindex $gene 0]
		}
		if {$len > [expr $mean+(8*$stdev)] || $len < [expr $mean-(8*$stdev)]} {
			lappend over_8stdev [lindex $gene 0]
		}
	}
	foreach gene $over_8stdev {
		db1 eval {SELECT binomial FROM t1 WHERE id = $gene} {
			set binomial "$binomial"
		}
		append over_8stdev_list "$binomial\t"
	}

	append output_line "\t$mean\t$stdev\t[llength $over_2stdev]\t[llength $over_8stdev]\t[string trim $over_8stdev_list]"
	append master_output "$output_line\n"
}

set out [open $direct/$org/$eval/Ortholog_fam_stats_$ival\.tsv w]
puts $out [string trim $master_output]
close $out
###########################################################################

db1 close
puts "DONE"




