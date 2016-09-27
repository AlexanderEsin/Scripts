### This proc takes the read input from a typical fasta files and splits it into a list of genes ###

proc split_genes {fasta} {
	global genes
	regsub -all {>} [string trim $fasta] {£>} fasta
	set fasta [split $fasta £]
	regsub -all "{}" $fasta {} genes
	return $genes
}