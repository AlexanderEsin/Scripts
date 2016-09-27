######
set direct ~/desktop/Geo_v_all/Geo_fastas
#####

proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}

###########################################################################

cd $direct

set geos [glob *.faa]

foreach geo $geos {
	openfile $geo
	set x [regexp -all {>} $data]
	set y [regsub -all {\(Bact\)} $data {(Geobac)} new_data]
	puts "$geo ==== $x ==== $y"
	if {$x != $y} {
		puts "ERROR"
	} else {
		set new_data [string trim $new_data]
		set out [open $geo w]
		puts $out $new_data
		close $out
	}
}