### CHECK WHETHER ANY OF THE GBFF DERIVED PROTEOMES CONTAIN A LOT OF MISSING PROTEIN IDS ###

source ~/Dropbox/Scripts/General_utils.tcl

cd /users/aesin/Desktop/Clean_Genomes/Gbff_derived_proteomes
set files [glob *faa]

foreach file $files {
	set data [string trim [openfile $file]]
	set genes [split_genes $data]
	set num_genes [llength $genes]

	set num_missing [regexp -all "missing_protein_id_qualifer" $data]
	spinner

	if {[expr $num_missing * 2] > $num_genes} {
		puts $file
	}

}

### RESULT: 

### ONLY 1 FILE: GCF_000234195.1_ASM23419v1_genomic.faa CONTAINS MORE THAN 50% MISSING IDS
### --> Candidatus Burkholderia kirkii UZHbot1 genomic scaffold