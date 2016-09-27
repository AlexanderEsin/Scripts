proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}

###########################################################################

cd ~/desktop/Geo_v_all/All_fasta_nogeobac

set files [glob *.faa]

foreach file $files {
	regsub -all " " $file {_} new_file
	file rename $file $new_file
}

set files [glob *.faa]

foreach file $files {
	openfile $file
	if {[regexp -all {=} $data] > 0} {
		set xx [regsub -all {=} $data {} new_data]
		set out [open $file w]
		puts $out $new_data
		close $out
		puts "$file ==== $xx substitutions"
	} else {
	}
}