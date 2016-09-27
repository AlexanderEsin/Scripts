#Proc to open DNA sequence fasta
proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}

#Open DNA sequence file. Split as before. Set i so that I can call the correct DNA sequence to
#match the alignment sequence by incrementing each cycle.

openfile Drosophila_notch.fa
			set nt_data $data
			regsub -all {>} $nt_data £> nt_data
			set glist [split $nt_data £]
			set glist [lrange $glist 1 end]
			set i 0

#Open sequence alignment file and split each gene.
openfile Muscle_Notch_aa.fa

			regsub -all {>} $data £> data
			set align_glist [split $data £]
			set align_glist [lrange $align_glist 1 end]

	     	foreach align_gene $align_glist {
				set aa_ori [string first \n $align_gene]
				set aa_seq [string range $align_gene $aa_ori end]

#4. Your script relies on having the genes and proteins in the same species order in the input file to muscle and the aligned output. Just be aware that the output need not be in the same order as the input! Depends on the program and the sequences. Muscle doesn't preserve order by default. So ideally, the script should contain some matching or prior sequence ordering.
				
				set align_header [string trim [string range $align_gene 0 $aa_ori]]

#Convert all the alignment dashes to X; clean out whitespace; append € to represent STOP

#5. converting dashes to X works here but it's superfluous really. Also the translate.tcl script converts stop codons to X and in some rare cases you might have a "stop" in the coding sequence (codes for selenocysteine in some proteins), in which case you'd run into trouble
	
				regsub -all {[-]} $aa_seq {X} aa_seq
				regsub -all {[^A-Z]} $aa_seq {} aa_seq
				append aa_seq {€}
			
#Set output file
				set dna ""

#Call the correct DNA sequence, then pull out the cds.

			
				set gene [lindex $glist $i]
				set ori [string first \n $gene]
				set cds [string range $gene $ori end]
				regsub -all {[^A-Z]} $cds {} cds

#Build output, misalignment noted as "X". € marker replaced by last three bases in each
#gene DNA sequence. Cds is only readusted if codon is used.

				while {[string length $aa_seq]>=1} {
					set aa [string index $aa_seq 0]
										
					if {[regexp {[X]} $aa]} {
						append dna {X}
					} elseif {[regexp {[€]} $aa]} {
						set codon [string range $cds end-2 end]
						append dna $codon
					} else {
						set codon [string range $cds 0 2]
						append dna $codon
						set cds [string range $cds 3 end]
					}
					set aa_seq [string range $aa_seq 1 end]
									}

#Alignment Xs replaced by *** to conserve frame and output looks cleaner than with dashes.
				incr i
				regsub -all {[^A-Z]} $dna {} dna
				regsub -all {[X]} $dna {---} dna
				set out [open Reverse_Transl.fa a]
				puts $out "$align_header\n$dna"
				close $out
			}
