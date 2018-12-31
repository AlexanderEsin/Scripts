#!/usr/local/bin/tclsh

source ~/Documents/Scripts/General_utils.tcl

set direct		/Users/aesin/Desktop/HTa_tonio
set align_dir	$direct/Alignment

set highHomo_seqs	$align_dir/highHomology.fasta
set lowHomo_seqs	$align_dir/lowHomology.fasta

# Align the high homology sequences in the first round
exec mafft-linsi $highHomo_seqs > $align_dir/highHomo_align.aln
exec mafft --add $lowHomo_seqs --reorder $align_dir/highHomo_align.aln > $align_dir/fullSeqs.aln


FastTreeMP -lg -gamma  $align_dir/fullSeqs.aln > fastTree_tree.txt

exec raxml -f a -s  $align_dir/fullSeqs.aln -n raxml_tree.txt -m PROTCATAUTO -p [expr {rand() * 10000}] -x [expr {rand() * 10000}] -N 100 -T 20