###
set direct ~/desktop/Geobac
###

proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}

###########################################################################

cd $direct

set prot_input ""
set patx {>+?.+?\]+?}

cd $direct/marked_faa
set protein_seqs [glob *.faa]
foreach species $protein_seqs {
	openfile $species
	set data [string trim $data]
	regsub -all {>} $data {£>} data
	set data [split $data £]
	regsub -all "{}" $data {} data
	foreach gene $data {
		regexp $patx $gene hit1
		append prot_input "$hit1\n"
		set out [open $direct/all_prots.txt a]
		puts $out $gene
		close $out
	}
}

set out2 [open $direct/faa_ids.txt a]
puts $out2 $prot_input
close $out2

###########################################################################


##This script takes all the individual fasta files (in faa_all) and produces two output files:
##faa_ids.txt contains just the headers of all the proteins in all the fasta files for faster label searching.
##all_prots.txt contains all the genes (with headers) for eventual sequence retrieval.