#!/bin/sh
#PBS -l walltime=25:00:00
#PBS -l mem=8gb
#PBS -l ncpus=16

cd $SCRATCH

echo "
set direct $SCRATCH/out_split

proc openfile {fl} {
	global data
	set in [open \$fl r]
	set data [read \$in]
	close \$in
	return
}

cd \$direct
set z [glob -nocomplain -type d *]
if {[regexp Rbbh \$z] == 1} {
} else {
	file mkdir Rbbh/Bin
}

cd \$direct/Out[lindex \$argv 0]
set blastlist [glob *.tsv]
set donelist \"\"
set RBBH_header \"GeneA\tGeneB\tEval(A-B)\tEval(B-A)\n\"
set hitnohit_header \"Organism_A\tOrganism_B\tTrue_Hit\tMiss_Hit\tProcessed\tNo_Blast_Hit(A)\tNo_Blast_Hit(B)\tTotal_Genes(A)\tTotal_Genes(B)\n\"
set hit_table \"\"
append hit_table \$hitnohit_header


foreach matchA \$blastlist {
	set rearrange [split \$matchA &]
	set q [lindex \$rearrange 0]
	set s [lindex \$rearrange 1]
	regsub .tsv \$s {} s
	set pat \$s&\$q
	append pat .tsv
	set matchlist [glob *.tsv]

	if {[regexp \$pat \$donelist] == 1} {
		puts \"Already Done!\"
	} elseif {[regexp \$pat \$matchlist matchB]} {
		openfile \$matchA
		set dataA \$data
		openfile \$matchB
		set dataB \$data
		append donelist \"\$matchA \"

		set q_nohit [regsub -all \"0 hits found\" \$dataA {} numb]
		set s_nohit [regsub -all \"0 hits found\" \$dataB {} numb]
		set tot_geneA [regsub -all \"hits found\" \$dataA {} numb]
		set tot_geneB [regsub -all \"hits found\" \$dataB {} numb]

		set besthit \"\"
		append besthit \$RBBH_header
		set hit_counter 0
		set miss_counter 0
		set processed 0
		set hitpat \"\$q\t\$s\"

		set lineparse [split \$dataA \"\n\"]
		regsub -all \"{}\" \$lineparse {} lineparse

		foreach line \$lineparse {
			if {[lindex \$line 0] == \"\#\"} {
			} else {
				set geneA [lindex \$line 0]
				set geneB [lindex \$line 1]
				set pata {\s+}
				set patb {\s+.+?}
				if {[regexp -line \$geneB\$pata\$geneA \$dataB match] == 1} {
					set evalA [lindex \$line 10]
					regexp -line \$geneB\$pata\$geneA\$patb \$dataB hitline
					set evalB [lindex \$hitline end-1]
					append besthit \"\$geneA\t\$geneB\t\$evalA\n\$geneB\t\$geneA\t\$evalB\n\"
					incr hit_counter
				} else {
				incr miss_counter
				}
			incr processed
			}
		}
		append hit_table \"\$hitpat\t\$hit_counter\t\$miss_counter\t\$processed\t\$q_nohit\t\$s_nohit\t\$tot_geneA\t\$tot_geneB\n\"
		set out [open \$direct/Rbbh/Bin/RBBH_\$q&\$s.txt a]
		puts \$out \"\$besthit\"
		close \$out
		puts \"Done:[llength \$donelist]\"
	} else {
		puts \"ERROR\"
	}
}
set out2 [open \$direct/Rbbh/RBBH_Hitfile.txt a]
puts \$out2 \"\$hit_table\"
close \$out2f
" > $SCRATCH/out_split/Rbbh_split_script.tcl

chmod +x $SCRATCH/out_split/Rbbh_split_script.tcl

tclsh8.5 $SCRATCH/out_split/Rbbh_split_script.tcl 1 &
tclsh8.5 $SCRATCH/out_split/Rbbh_split_script.tcl 2 &
tclsh8.5 $SCRATCH/out_split/Rbbh_split_script.tcl 3 &
tclsh8.5 $SCRATCH/out_split/Rbbh_split_script.tcl 4 &
tclsh8.5 $SCRATCH/out_split/Rbbh_split_script.tcl 5 &
tclsh8.5 $SCRATCH/out_split/Rbbh_split_script.tcl 6 &
tclsh8.5 $SCRATCH/out_split/Rbbh_split_script.tcl 7 &
tclsh8.5 $SCRATCH/out_split/Rbbh_split_script.tcl 8 &
tclsh8.5 $SCRATCH/out_split/Rbbh_split_script.tcl 9 &
tclsh8.5 $SCRATCH/out_split/Rbbh_split_script.tcl 10 &
tclsh8.5 $SCRATCH/out_split/Rbbh_split_script.tcl 11 &
tclsh8.5 $SCRATCH/out_split/Rbbh_split_script.tcl 12 &
tclsh8.5 $SCRATCH/out_split/Rbbh_split_script.tcl 13 &
tclsh8.5 $SCRATCH/out_split/Rbbh_split_script.tcl 14 &
tclsh8.5 $SCRATCH/out_split/Rbbh_split_script.tcl 15 &
tclsh8.5 $SCRATCH/out_split/Rbbh_split_script.tcl 16 &

wait

cd $SCRATCH
rm *tcl


