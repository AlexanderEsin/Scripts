set direct /users/aesin/desktop/Clean_proteomes

set group All_fastas_geo_mark

set output_dir /users/aesin/desktop/Clean_proteomes

###########################################################################

cd $direct/$group
set fastas [glob *faa]
set i 1

foreach fasta $fastas {
	openfile $fasta
	set data [string trim $data]
	set out [open $output_dir/$group\_complete.faa a]
	puts $out $data
	close $out

	puts "Combining fasta files ... $i/[llength $fastas]"
	incr i
}
