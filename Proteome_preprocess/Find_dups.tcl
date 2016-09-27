#!/usr/local/bin/tclsh
###
set org Viral
set direct /users/aesin/desktop/Cleaning_prots/$org
###

proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}

# This proc remakes the 80-wide fasta entries from a flat CDS sequence #
proc linebreak {s {width 80}} {
   global res
   set res ""
   while {[string length $s]>$width} {
       set res "$res[string range $s 0 79]\n"
       set s [string range $s 80 end]
   }
   set res "$res$s"
   return $res
}

proc split_genes {fasta} {
	global genes
	regsub -all {>} [string trim $fasta] {£>} fasta
	set fasta [split $fasta £]
	regsub -all "{}" $fasta {} genes
}

###
# Make the necessary subfolders #
cd $direct
file mkdir MSA_check
file mkdir Duplicates
file mkdir Clean_faa

###########################################################################

cd $direct/Orig_faa
# List all the present fasta files in the directory #
set fastas [glob *.faa]
set fasta_no [llength $fastas]
set y 1

set fasta [lindex $fastas 0]
foreach fasta $fastas  {
	set flat_glist {}
	set dupe_id_list {}
	set duplicate_counter 1

	openfile $fasta

	# Make a list of genes for each fasta, and for each gene flatten the CDS (take out the newlines) - creating a new gene list #
	split_genes [string trim $data]
	set gene_no [llength $genes]
	foreach gene $genes {
		set header [string range $gene 0 [string first \n $gene]]
		set CDS [string range $gene [string first \n $gene] end]
		regsub -all {\n} $CDS {} flat_CDS
		lappend flat_glist "$header$flat_CDS"
	}

	# For each gene in the new, flattened gene list we search against the rest of the genes for matches #
	while {[llength $flat_glist] > 0} {

		set test_gene [lindex $flat_glist 0]
		set test_CDS [lindex [split $test_gene \n] 1]
		# A list of all the hits #
		set hit_indices [lsearch -all -glob $flat_glist *\n$test_CDS]

		# If no hits are found, that means the CDS did not find itself, which is an error #
		if {[llength $hit_indices] == 1 && $hit_indices == -1} {
			puts "ERROR: the searched test_CDS could not find itself..."
			puts "\nFile: $fasta\nGene: $test_gene\n"
			exit 2
		} else {
			# Copy number should be at least 1 (found itself) #
			set copy_number [llength $hit_indices]
		}

		# If duplicates were found... #
		if {$copy_number > 1} {
			# Make a new folder to hold the duplicate sequences #
			file mkdir $direct/MSA_check/$fasta/$duplicate_counter
			set duplicate_list {}
			set temp_header_list {}
			# For each duplicate hit (including the original query) #
			foreach index $hit_indices {
				set duplicate [lindex $flat_glist $index]
				set header [string trim [string range $duplicate 0 [string first \n $duplicate]]]
				set flat_CDS [string trim [string range $duplicate [string first \n $duplicate] end]]

				set reconst_CDS [linebreak $flat_CDS]
				set reconst_dup_gene "$header\n$reconst_CDS"

				lappend duplicate_list $reconst_dup_gene
				lappend temp_header_list $header
			}
			set out [open $direct/MSA_check/$fasta/$duplicate_counter/$duplicate_counter\.faa w]
			puts $out [string trim [join $duplicate_list \n]]
			close $out
			incr duplicate_counter
			lappend dupe_id_list [join $temp_header_list \n]
		}

		if {$copy_number > 2} {
			set extra_hits [lrange $hit_indices 1 end-1]
			set adjustment 0
			foreach extra_hit $extra_hits {
				set flat_glist [lreplace $flat_glist [expr $extra_hit - $adjustment] [expr $extra_hit - $adjustment]]
				incr adjustment
			}
		}
		puts "$y/$fasta_no == [expr $gene_no -[llength $flat_glist]]/$gene_no"
		set flat_glist [lrange $flat_glist 1 end]
	}
	if {[llength $dupe_id_list] > 0} {
		set out [open $direct/Duplicates/$fasta.txt w]
		puts $out [join $dupe_id_list "\n\n"]
		close $out
	}
	incr y
}

