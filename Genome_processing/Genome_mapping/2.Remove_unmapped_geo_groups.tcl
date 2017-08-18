#!/usr/local/bin/tclsh
source ~/Dropbox/Scripts/General_utils.tcl
package require sqlite3

## Set paths ##
set direct /users/aesin/desktop
set working_directory $direct/Geo_analysis/Geo_ortholog_nucl

set output_directory $working_directory/Groups_all_geo_map
set output_fail_directory $working_directory/Groups_not_all_geo_map

file mkdir $output_directory
file mkdir $output_fail_directory

## Databases ##
sqlite3 all_prot_geo_db $direct/Clean_Proteomes/all_prot_geo_ncbi_db
sqlite3 genome_posit_db $working_directory/Geo_v_all_orth_location_db

cd $working_directory/Groups

set fasta_files [lsort -dictionary [glob *faa]]
progress_init [llength $fasta_files]

set total_geobac 0
set total_present 0
set i 0

###########################################################################

foreach file $fasta_files {
	set prot_seqs [split_genes [openfile $file]]
	set prot_seqs_original $prot_seqs
	set group [string range $file 0 end-4]

	set prot_to_remove {}

	## Identify which Geobacillus proteins do not have mappable coordinates ##
	foreach prot_seq $prot_seqs {
		set prot_seq [string trim $prot_seq]
		set prot_id [string range $prot_seq 1 [string first " " $prot_seq]-1]
		
		set class [string range [all_prot_geo_db eval {SELECT class FROM t1 where id = $prot_id}] 1 end-1]
		
		if {$class eq "Geobac"} {
			set present [genome_posit_db eval {SELECT EXISTS(SELECT 1 FROM t1 WHERE prot_id = $prot_id LIMIT 1)}]

			set total_present [expr $total_present + $present]
			incr total_geobac

			if {$present == 0} {
				lappend prot_to_remove $prot_id
			} else {
				continue
			}
		}
	}

	## Remove these from the orthologous protein group ##
	foreach protein $prot_to_remove {
		set idx_list [lsearch -all -glob $prot_seqs *$protein*]
		foreach idx [lreverse $idx_list] {
			set prot_seqs [lreplace $prot_seqs $idx $idx]
		}
	}

	## Check whether the family is still eligible for inclusion ##
	# Are there still four or more proteins? #
	if {[llength $prot_seqs] < 4} {
		puts "Group: $group no longer has at least 4 proteins\n"
		lappend removed_groups $group
		set out [open $output_fail_directory/$file w]
		puts $out [join $prot_seqs_original ""]
		close $out
	} elseif {[llength [lsearch -all -glob $prot_seqs *\(Geobac\)*]] == 0} {
		puts "Group: $group no longer has at least one Geobacillus protein\n"
		lappend removed_groups $group
		set out [open $output_fail_directory/$file w]
		puts $out [join $prot_seqs_original ""]
		close $out
	} else {
		set out [open $output_directory/$file w]
		puts $out [join $prot_seqs ""]
		close $out
	}

	incr i
	progress_tick $i
}

puts "Total Geobacilli queried: $total_geobac"
puts "Total with coordinates: $total_present"
puts "Total groups removed: [llength $removed_groups]"

set out [open $working_directory/Removed_groups.tsv w]
puts $out [join $removed_groups \n]
close $out



#puts "\n\n\n\n$total_geobac"
all_prot_geo_db close
genome_posit_db close
