

set direct /users/aesin/desktop/Geo_v_all/2.0
set group_type Anoxy_geo_only_groups

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

cd $direct/$group_type/Group_fastas

set group_files [glob *faa]
set output_table {}


foreach file $group_files {

	set total_aa_count 0

	openfile $file
	split_genes $data

	foreach gene $genes {
		set CDS [string trim [string range $gene [string first \n $gene] end]]
		regsub -all \n $CDS {} flat_CDS
		set flat_CDS [split [string trim $flat_CDS] ""]
		set aa_count [llength $flat_CDS]
		set total_aa_count [expr $total_aa_count + $aa_count]
	}

	set average_aa_count [expr round($total_aa_count / double([llength $genes]))]
	puts "File: $file\tAverage aa count: $average_aa_count"
	lappend output_table $average_aa_count
}

set out [open $direct/$group_type/Average_aa_counts.txt w]
puts $out [join $output_table \n]
close $out