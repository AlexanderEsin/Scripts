#!/usr/local/bin/tclsh
###########################################################################
### Changelog ###
# 14 January 2016 - Complete revamp. Cut 50+ lines of unnecessary code, converted many variables to lists for speed increase #


###########################################################################
### Procedures and packages ###

source ~/Dropbox/Scripts/General_utils.tcl
package require sqlite3
###########################################################################
### Set variables ###

set ival 2.0
set direct /users/aesin/desktop/Geo_analysis/Geo_v_all/$ival; # For doing the missing Geobacillus families: /users/aesin/desktop/Geo_analysis/Geo_v_all/$ival/Missing_geobac_reconstruction #
# set direct /users/aesin/desktop/Geo_analysis/Geo_v_all/$ival/Missing_geobac_reconstruction
set input_fasta_directory Fastas

set split_core_switch 1
set perc_cutoff 25

###########################################################################
### Load the database ###
sqlite3 db1 /users/aesin/desktop/Clean_proteomes/all_prot_geo_ncbi_db

###########################################################################
### If there is at least a single Geobacillus per group, construct a core group based on the Geobacillus gene length with best coverage. Output to Core_Aligns ###

file mkdir $direct/Group_fastas_MSA
file mkdir $direct/Core_aligns

cd $direct/$input_fasta_directory

set y 1
set geo_core_comparison [list "Group\tGroup\tOrig_prots\tCore_prots\tOriginal_geo_prots\tCore_geo_prots\tCore_CDS_length\tCore_CDS_ID"]
set fastas [lsort -dictionary [glob *.faa]]
set perc_cutoff [expr double($perc_cutoff) / 100]

foreach fasta $fastas {

	set group [string range $fasta 0 end-4]
	file mkdir $direct/Core_aligns/$group 

	set measures {}

	####################################
	set fasta_data [openfile $fasta]
	set genes [split_genes $fasta_data]

	set geobac_headers [regexp -all -inline -line {.+\(Geobac\)+.+} $fasta_data]
	set orig_geo_no [llength $geobac_headers] 

	foreach header $geobac_headers {
		set prot_id [string range $header 1 [string first " " $header]-1]
		
		## In this request the output is a list of 1, so we select the first index ##
		set gene [lindex [db1 eval {SELECT gene FROM t1 WHERE id = $prot_id}] 0]

		regsub -all "\n" [string range $gene [string first \n $gene]+1 end] {} CDS	
		lappend measures [list $prot_id [string trim $CDS]]
	}

	####################################
	### If the number of Geobacillus proteins is 0 - do not filter into a core genome, and produce a key for the genes ###
	if {[llength $measures] == 0} {
		set i 1
		set new_glist {}
		set key_list {}

		foreach gene $genes {
			## Get the protein ID, class, and binomial of each gene ##
			set id [string range $gene 1 [string first " " $gene]-1]
			set class [string range $gene [string first \( $gene]+1 [string first \) $gene]-1]
			set binomial [string range $gene [string last \[ $gene]+1 [string first \n $gene]-2]

			## Replace the header with the k(i) nomenclature for alignment ##
			regsub {>+?.+?\n+?} $gene ">k$i\n" new_gene
			lappend new_glist [string trim $new_gene]
			lappend key_list "k$i\t$binomial\t$class\t$id"
			incr i
		}

		set out [open $direct/Core_aligns/$group/CORE_$fasta w]
		puts $out [join $new_glist \n]
		close $out

		set out [open $direct/Core_aligns/$group/KEY_$group.txt w]
		puts $out [join $key_list \n]
		close $out

		####################################

		lappend geo_core_comparison "$group\tORIGINAL\t$orig_prot_no\t$orig_prot_no\t$orig_geo_no\t$orig_geo_no\t-\t-"

	### Otherwise, for each geobacillus CDS, see how many of the other CDS will be included given +/- $perc_cutoff length boundaries on that CDS. Pick the CDS that includes the most. ###
	} else {
		####################################
		set win_no 0

		foreach query_measure $measures {
			set strlength [string length [lindex $query_measure 1]]
			set high_cut_off [::tcl::mathfunc::round [expr $strlength+($strlength * $perc_cutoff)]]
			set low_cut_off [::tcl::mathfunc::round [expr $strlength-($strlength * $perc_cutoff)]]

			## Count how many Geobacilli would be included in +/- 25% of sequence length of each Geobacillus ##
			set included_geo 0
			foreach subject_measure $measures {
				if {[string length [lindex $subject_measure 1]] <= $high_cut_off && [string length [lindex $subject_measure 1]] >= $low_cut_off} {
					incr included_geo
				}
			}

			## If the number of other Geobacilli included within the boundary is greater than for the previous best, set current as the best ##
			if {$included_geo > $win_no} {
				set win_no $included_geo
				set best_strlength $strlength
				set win_measure [lindex $query_measure 0]

			## If the number of Geobacilli included is equal (e.g. there are only 2, and they do not fit into each others' brackets) then pick the one with the shortest sequence ##
			} elseif {$included_geo == $win_no} {
				if {$strlength < $best_strlength} {
					set win_no $included_geo
					set best_strlength $strlength
					set win_measure [lindex $query_measure 0]
				}
			}
		}

		####################################
		### Set some variables for downstream ###

		### Get the identity (gene) of the winner for sorting downstream ###
		set win_id $win_measure

		### Set the +/- $perc_cutoff length boundaries on the best CDS as calculated above ###
		set high_cut_off [::tcl::mathfunc::round [expr $best_strlength+($best_strlength * $perc_cutoff)]]
		set low_cut_off [::tcl::mathfunc::round [expr $best_strlength-($best_strlength * $perc_cutoff)]]

		### Filter the rest of the group for proteins that lie within the set boundaries to get the "core" genome, and produce a key that can be used for the entire group ###
		set core_counter 0
		set other_counter 0

		set core_fasta {}
		set other_fasta {}

		set new_glist {}
		set key_list {}
		set geobac_k_list {}

		####################################
		### Produce a key - it will be the same regardless of sequence inclusion in CORE or OTHER sets ###
		set i 1
		foreach gene $genes {
			## Get the protein ID, class, and binomial of each gene ##
			set id [string range $gene 1 [string first " " $gene]-1]
			set class [string range $gene [string first \( $gene]+1 [string first \) $gene]-1]
			set binomial [string range $gene [string last \[ $gene]+1 [string first \n $gene]-2]

			## Replace the header with the k(i) nomenclature for alignment ##
			regsub {>+?.+?\n+?} $gene ">k$i\n" new_gene
			lappend new_glist [string trim $new_gene]
			lappend key_list "k$i\t$binomial\t$class\t$id"
			
			## We wil count how many Geobacillus are included in the core ##
			if {[string first "Geobacillus" $gene] > 0} {
				lappend geobac_k_list ">k$i"
			}

			incr i
		}
		
		set out [open $direct/Core_aligns/$group/KEY_$group.txt w]
		puts $out [join $key_list \n]
		close $out

		####################################

		# set orig_prot_no [regsub -all {>} $new_glist {£>} genes]
		# set genes [split [string trim $genes] £]
		# set genes [lrange $genes 1 end]
		set orig_prot_no [llength $genes]
		set genes $new_glist

		### Filter the protein sequences into either the "core" or "other" sets depending on whether they are within +/-25% of the sequence length of the representative Geobacillus sequence ###
		foreach gene $genes {
			regsub -all "\n" [string range $gene [string first \n $gene]+1 end] {} CDS	
			if {[string length $CDS] <= $high_cut_off && [string length $CDS] >= $low_cut_off} {
				# append core_fasta "$gene"
				lappend core_fasta $gene
				incr core_counter
			} else {
				# append other_fasta "$gene"
				lappend other_fasta $gene
				incr other_counter
			}
		}

		### If the number of proteins pulled into the core genome is 1 (i.e only the Geobacillus gene), then just use the original fasta ###
		if {$core_counter == 1} {
			set out [open $direct/Core_aligns/$group/CORE_$fasta w]
			puts $out [join $new_glist \n]
			close $out

			lappend geo_core_comparison "$group\tORIGINAL\t$orig_prot_no\t$orig_prot_no\t$orig_geo_no\t$orig_geo_no\t$best_strlength\t$win_id\n"

		### Otherwise count the number of Geobacilli included in the new core genome and output it ###
		} else {
			set out [open $direct/Core_aligns/$group/CORE_$fasta w]
			puts $out [join $core_fasta \n]
			close $out

			## Count the number of Geobacillus that made it into the core. This should of course be the same number as the winning "geo_included" from above ##
			set core_geo_no 0
			foreach geo_k_val $geobac_k_list {
				if {[lsearch -glob $core_fasta $geo_k_val\n*] != -1} {
					incr core_geo_no
				}
			}
			
			if {$other_counter > 0} {
				set out [open $direct/Core_aligns/$group/OTHER_$fasta w]
				puts $out [join $other_fasta \n]
				close $out

				### Also output an win_id_X file that contains the id of the core geobac sequence and its CDS length
				set out [open $direct/Core_aligns/$group/win_id_$group w]
				puts $out "Win_id:\t[string trim $win_id]\nWin_length:\t$best_strlength"
				close $out
			}

			lappend geo_core_comparison "$group\tCORE\t$orig_prot_no\t$core_counter\t$orig_geo_no\t$core_geo_no\t$best_strlength\t$win_id\n"
		}		
	}
	puts "Done: $y / [llength $fastas]"
	incr y
}

set out [open $direct/Geo_core_comparison.tsv w]
puts $out [join $geo_core_comparison \n]
close $out

db1 close

##########################################################################
### Split the files into n directories and then run MSA on directories (see MSA_master and MSA_slave scripts for details). Clustal Omega runs on 1 core per alignment. ###

if {$split_core_switch == 1} {
	cd $direct
	file mkdir Core_aligns_split
	#n=number of cores wanted - 1 (e.g. for 64 cores n = 63; max 7 (i.e. 8) for quad-core desktop MAC)
	set n 59
	###

	cd $direct/Core_aligns
	set org_size_list [lsort -dictionary [exec ls]]

	set i 0
	set y 1

	while {$y <= [llength $org_size_list]} {
		while {$i <= $n} {
			set dir [expr ($i + 1)]
			file mkdir $direct/Core_aligns_split/$dir
			set counter $i
			while {$counter < [llength $org_size_list]} {
				file copy -force [lindex $org_size_list $counter] $direct/Core_aligns_split/$dir/
				set counter [expr ($counter + ($n + 1))]
				incr y
			}
			incr i
		}
	}
}

##########################################################################