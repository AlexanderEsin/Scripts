source ~/Dropbox/Scripts/General_utils.tcl
package require sqlite3

###########################################################################

cd /users/aesin/desktop/Clean_proteomes
sqlite3 db1 all_prot_geo_db

###########################################################################

set direct /users/aesin/Desktop/Geo_analysis/Geo_omes
cd $direct/Geobac_proteomes

###########################################################################

set geobac_proteomes [glob *faa]
set tabular_out {}
lappend tabular_out "Name\tNumber of genes"
foreach proteome_file $geobac_proteomes {

	set genes_fasta [string trim [openfile $proteome_file]]
	set genes [split_genes $genes_fasta]

	set sample_gene [lindex $genes 1]
	set sample_id [string range $sample_gene 1 [string first " " $sample_gene]-1]

	db1 eval {SELECT binomial FROM t1 WHERE id = $sample_id} {
		set binomial "$binomial"
	}

	lappend tabular_out "$binomial\t[llength $genes]"
}

set out [open $direct/Geobac_num_proteins.tsv w]
puts $out [join $tabular_out \n]
close $out


db1 close