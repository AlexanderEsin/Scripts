set direct /users/aesin/desktop
set group_type Nevena_phylo

## PROCS ##
### Open and return the contents of a file ###
proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}

### Split a fasta file full of proteins into a list of individual proteins ###
proc split_genes {fasta} {
	global genes
	regsub -all {>} [string trim $fasta] {£>} fasta
	set genes [lrange [split $fasta £] 1 end]
    return
}

proc tcl::mathfunc::roundto {value decimalplaces} {
	expr {round(10**$decimalplaces*$value)/10.0**$decimalplaces}
}

###########################################################################

cd $direct/$group_type

set fastas [glob *txt]

set total_genome {}

foreach fasta $fastas {
	openfile $fasta
	split_genes $data
	foreach gene $genes {
		lappend total_genome $gene
	}
}


###########################################################################

puts "Total number of genes: [llength $total_genome]"

set aa_counts {}
set i 1

foreach gene $total_genome {
	set CDS [string trim [string range $gene [string first \n $gene] end]]
	regsub -all \n $CDS {} flat_CDS
	set flat_CDS [split [string trim $flat_CDS] ""]
	set aa_count [llength $flat_CDS]
	lappend aa_counts $aa_count
	incr i
}

set aa_counts [lsort -increasing -dictionary $aa_counts]

set out [open $direct/$group_type/Total_aa_counts.tsv w]
puts $out [join $aa_counts \n]
close $out
