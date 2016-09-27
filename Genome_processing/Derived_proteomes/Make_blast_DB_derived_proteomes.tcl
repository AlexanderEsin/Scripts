#!/usr/local/bin/tclsh
set direct /users/aesin/desktop/Clean_Genomes
file mkdir $direct/Gbff_derived_proteomes_db

cd $direct/Gbff_derived_proteomes
set files [glob *faa]
set num_files [llength $files]
set i 0

foreach file $files {
	incr i
	set trimmed_name [string range $file 0 end-4]
	set db_name "$trimmed_name\_db"

	catch {exec makeblastdb -in $file -dbtype prot -out ../Gbff_derived_proteomes_db/$db_name}
	puts "Done: $i / $num_files"
}