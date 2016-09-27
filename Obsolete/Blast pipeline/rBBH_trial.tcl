#Change the directory below
set direct ~/Desktop/Fungi-Archaea-Bact/Klast_Bacteria_All

proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}

cd $direct
set z [glob -nocomplain -type d *]
if {[regexp Rbbh $z] == 1} {
} else {
	file mkdir Rbbh/Bin
}

#List a set of .tsv files. Donelist used to avoid reverse rbbh-ing.

cd $direct/out
set blastlist [glob *.tsv]
set donelist ""
set RBBH_header "GeneA\tGeneB\tEval(A-B)\tEval(B-A)\n"
set hitnohit_header "Organism_A\tOrganism_B\tTrue_Hit\tMiss_Hit\tProcessed\tNo_Blast_Hit(A)\tNo_Blast_Hit(B)\tTotal_Genes(A)\tTotal_Genes(B)\n"
set hit_table ""
append hit_table $hitnohit_header

#For each file in $blastlist, rearrange the template to match the reciprocal file. Matchlist contains the same .tsv files for comparison with matchA.

foreach matchA $blastlist {
	set rearrange [split $matchA &]
	set q [lindex $rearrange 0]
	set s [lindex $rearrange 1]
	regsub .tsv $s {} s
	set pat $s&$q
	append pat .tsv
	set matchlist [glob *.tsv]

	#As long as the reciprocal hasn't been done yet (regexp below), file_A and reciprocal are opened, and added to $donelist.

	if {[regexp $pat $donelist] == 1} {
		puts "Already Done!"
	} elseif {[regexp $pat $matchlist matchB]} {
		openfile $matchA
		set dataA $data
		openfile $matchB
		set dataB $data
		append donelist "$matchA "

		#0 hits and total number of genes metrics are defined here for each file

		set q_nohit [regsub -all "0 hits found" $dataA {} numb]
		set s_nohit [regsub -all "0 hits found" $dataB {} numb]
		set tot_geneA [regsub -all "hits found" $dataA {} numb]
		set tot_geneB [regsub -all "hits found" $dataB {} numb]

		#Besthit will hold output.

		set besthit ""
		append besthit $RBBH_header
		set hit_counter 0
		set miss_counter 0
		set processed 0
		set hitpat "$q\t$s"

		set lineparse [split $dataA "\n"]
		regsub -all "{}" $lineparse {} lineparse

		foreach line $lineparse {
			if {[lindex $line 0] == "\#"} {
			} else {
				set geneA [lindex $line 0]
				set geneB [lindex $line 1]
				set pata {\s+}
				set patb {\s+.+?}
				if {[regexp -line $geneB$pata$geneA $dataB match] == 1} {
					set evalA [lindex $line 10]
					regexp -line $geneB$pata$geneA$patb $dataB hitline
					set evalB [lindex $hitline end-1]
					append besthit "$geneA\t$geneB\t$evalA\n$geneB\t$geneA\t$evalB\n"
					incr hit_counter
				} else {
				incr miss_counter
				}
			incr processed
			}
		}
		append hit_table "$hitpat\t$hit_counter\t$miss_counter\t$processed\t$q_nohit\t$s_nohit\t$tot_geneA\t$tot_geneB\n"
		set out [open $direct/Rbbh/Bin/RBBH_$q&$s.txt a]
		puts $out "$besthit"
		close $out
		puts "Done:[llength $donelist]"
	}
}
set out2 [open $direct/Rbbh/RBBH_Hitfile.txt a]
puts $out2 "$hit_table"
close $out2
