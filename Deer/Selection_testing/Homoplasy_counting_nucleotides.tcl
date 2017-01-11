## Homoplast detection - WIP ##

#!/usr/local/bin/tclsh

source ~/Dropbox/Scripts/General_utils.tcl

cd /users/aesin/Desktop/Deer/Selection_analysis/Alignments_new/Aligned

set align_data [string trim [openfile Aligned_HBB_names_only_ruminants.fas]]
set genes [split_genes_fast $align_data]

proc split_by_codon {gene} {
	set codon_l {}
	while {[string length $gene] > 3} {
		set codon [string range $gene 0 2]
		lappend codon_l $codon
		set gene [string range $gene 3 end]
	}
	return $codon_l
}

proc split_codon_by_nucl {codon_l} {

	set pos_1 {}
	set pos_2 {}
	set pos_3 {}

	foreach codon $codon_l {
		set nucls [split $codon ""]
		lappend pos_1 [lindex $nucls 0]; lappend pos_2 [lindex $nucls 1]; lappend pos_3 [lindex $nucls 2]
	}
	set pos_1 [join [lsort -unique $pos_1] " "]; set pos_2 [join [lsort -unique $pos_2] " "]; set pos_3 [join [lsort -unique $pos_3] " "]
	set var_by_pos "$pos_1\n$pos_2\n$pos_3"
	return $var_by_pos
}


## Process the alignment into codons ##
set alignment_by_codon {}
foreach gene $genes {
	set gene_by_codon_l {}
	set CDS [string trim [string range $gene [string first \n $gene] end]]
	regsub -all "\n" $CDS {} CDS_no_gap

	set gene_by_codon [split_by_codon $CDS_no_gap]
	
	for {set i 1} {$i <= [llength $gene_by_codon]} {incr i} {
		set codon [lindex $gene_by_codon [expr $i - 1]]

		if {$gene eq [lindex $genes 0]} {
			set codon_$i {}
		}

		lappend codon_$i $codon

		if {$gene eq [lindex $genes end]} {
			lappend alignment_by_codon [join [expr $\codon_$i] \t]
		}
	}
}

set codon_counter 1
foreach codons $alignment_by_codon {
	set variation [lsort -unique $codons]

	if {[llength $variation] == 1} {
		#puts "No variation at site $codon_counter"
	} else {
		set nuc_var [split_codon_by_nucl $codons]
		#puts "Site $codon_counter:\n$nuc_var\n\n"
		set by_pos [split $nuc_var \n]

		set var_in_sites 0
		for {set i 0} {$i <= [llength $by_pos]} {incr i} {
			set pos [lindex $by_pos $i]
			set variation [llength [split [string trim $pos] " "]]
			if {$variation > 1} {
				incr var_in_sites
			}
		}
		if {$var_in_sites > 2} {
			puts "Variation at $var_in_sites positions at site $codon_counter"
		}
		
	}
	incr codon_counter
}

