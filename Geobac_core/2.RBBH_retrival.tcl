# Change the directories below #
set direct /users/aesin/Desktop/Consensus_trees/Geobac
set org All_geobac

# A cut off value for which results to pick from the blast output - the default blast has a cut off of 1e-10 #
set cut_off 1e-100
regsub {1e} $cut_off {} folder_eval

proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}

###########################################################################

# Make output folder #
cd $direct/$org
file mkdir Rbbh_all_$folder_eval/Bin

# Prepare output variables #
set RBBH_header "GeneA\tGeneB\tEval(A-B)\tEval(B-A)\n"
set hitnohit_header "Organism_A\tOrganism_B\tTrue_Hit\tBelow_cutoff\tMiss_Hit\tProcessed\tNo_Blast_Hit(A)\tNo_Blast_Hit(B)\tTotal_Genes(A)\tTotal_Genes(B)\n"
set hit_table ""

# Go to the relevant blast hit directory and set regex patterns #
cd $direct/$org/Out_all
set orig_list [glob *.tsv]
set blastlist $orig_list
file mkdir temp

set pata {\s+}
set patb {\s+.+?}
set patc {\#+?.+?\n+?}

set i 1

# For each blast hit, find the reciprocal & open both #
while {[llength $blastlist] > 0} {
	set matchA [lindex $blastlist 0]
	set rearrange [split $matchA &]
	set q [lindex $rearrange 0]
	set s [lindex $rearrange 1]
	regsub .tsv $s {} s
	
	set matchB "$s&$q\.tsv"
	openfile $matchA
	set dataA $data
	openfile $matchB
	set dataB $data

	# Count the total gene number and no-hits for both files. The "n" subVar is necessary, as without it the regsub changes function #
	set q_nohit [regsub -all "0 hits found" $dataA {} n]
	set s_nohit [regsub -all "0 hits found" $dataB {} n]
	set tot_geneA [regsub -all "hits found" $dataA {} n]
	set tot_geneB [regsub -all "hits found" $dataB {} n]

	# Besthit will hold output + set counters #
	set besthit "$RBBH_header"
	set hit_counter 0
	set below_cutoff 0
	set miss_counter 0
	set processed 0
	set hitpat "$q\t$s"

	# Remove all comment lines #
	regsub -all $patc $dataA {} dataA
	regsub -all $patc $dataB {} dataB
	set lineparse [string trim [split [string trim $dataA] \n]]
	set dataB [string trim $dataB]
	# For every hit in fileA, find the reciprocal (if true best hit). Pull out evalues and append the output variable with data. #
	foreach line $lineparse {
		set geneA [lindex $line 0]
		set geneB [lindex $line 1]
		set pos_match [regexp -line -inline $geneB$pata$geneA$patb $dataB]
		if {[llength $pos_match] == 1} {
			set pos_match [lindex $pos_match 0]
			set evalA [lindex $line 10]
			set evalB [lindex $pos_match end-1]
			# Screen the blast output to only capture results below a certain e-value threshold - a value of 1e-10 is the default (same as the blast output) #
			if {$cut_off == 1e-10} {
				append besthit "$geneA\t$geneB\t$evalA\n$geneB\t$geneA\t$evalB\n"
				incr hit_counter
			} else {
				if {$evalA < $cut_off && $evalB < $cut_off} {
					append besthit "$geneA\t$geneB\t$evalA\n$geneB\t$geneA\t$evalB\n"
					incr hit_counter
				} else {
					incr below_cutoff
				}
			}
		} else {
		incr miss_counter
		}
	incr processed
	}

	# Create output file #
	append hit_table "$hitpat\t$hit_counter\t$below_cutoff\t$miss_counter\t$processed\t$q_nohit\t$s_nohit\t$tot_geneA\t$tot_geneB\n"
	set out [open $direct/$org/Rbbh_all_$folder_eval/Bin/RBBH_$q&$s.txt w]
	puts $out [string trim $besthit]
	close $out
	
	# Remove both files (temporarily) from the list to increase search speeds for matchB (above) #
	file rename $matchA temp/$matchA
	file rename $matchB temp/$matchB

	puts "Done: $i / [expr ([llength $orig_list]/2)]"

	# Reappraise blastlist #
	set blastlist [glob -nocomplain *.tsv]
	incr i
}
# Print stats to global file, if the file does not yet exist - add a header #
cd $direct/$org/Rbbh_all_$folder_eval
if {[file exists RBBH_Hitfile.txt] == 0} {
	set hit_table "$hitnohit_header$hit_table"
}
set out2 [open $direct/$org/Rbbh_all_$folder_eval/RBBH_Hitfile_$folder_eval.txt a]
puts $out2 [string trim $hit_table]
close $out2
cd $direct/$org/Out_all

# Refill the original out directory with file from temp #
cd temp
set files [glob *tsv]
foreach file $files {
	file rename $file $direct/$org/Out_all/$file
}
cd ..
file delete -force temp

clear
puts "RBBH Retrival == DONE"



