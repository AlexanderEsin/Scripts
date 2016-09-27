# Script to label each gene_id in each proteome with a species_specfic number (found in key_file). This allows the database to be redundant and means that the same ID will no longer link to multiple species/proteomes #

set direct /users/aesin/desktop/Clean_proteomes
set folder Viral

proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}

proc split_genes {fasta} {
	global genes
	regsub -all {>} [string trim $fasta] {£>} fasta
	set fasta [split $fasta £]
	regsub -all "{}" $fasta {} genes
}

###########################################################################

cd $direct/$folder
set fastas [glob *faa]

set key_file "File Name\tID_suffix\n"

set i 1
foreach fasta $fastas {
	openfile $fasta
	split_genes $data

	set new_glist {}

	foreach gene $genes {
		set id [string trim [string range $gene 1 [string first " " $gene]]]
		set new_id "$id\.$i"
		regsub $id $gene $new_id new_gene
		lappend new_glist $new_gene
	}

	set new_glist [join $new_glist {}]

	set out [open $fasta w]
	puts $out [string trim $new_glist]
	close $out

	append key_file "$fasta\t$i\n"

	puts "$i / [llength $fastas]"
	incr i
}

set out [open $direct/$folder\_key_file.tsv w]
puts $out [string trim $key_file]
close $out