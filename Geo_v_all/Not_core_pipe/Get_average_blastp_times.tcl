set direct /users/aesin/desktop/Geo_v_all/Timing

proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}

###########################################################################

cd $direct

set time_files [glob *txt]
foreach file $time_files {
	openfile $file
	set new_time_out ""
	regsub -all {\n\n} $data "\n" data
	set times [split [string trim $data] \n]
	foreach time $times {
		set time [split $time \t]
		set time_1 [lindex $time 1]
		set time_2 [lindex $time 2]

		if {$time_1 < 0} {
			puts "$file $time_1"
			set time_1 2000
			puts "$file $time_1"
		}
		if {$time_2 < 0} {
			puts "$file $time_2"
			set time_2 2000
			puts "$file $time_2"
		}
		set average_t [expr ($time_1+$time_2)/2]
		append new_time_out "[lindex $time 0]\t$average_t\n"
	}
	regsub {.faa} $file {} file
	set out [open Average_$file w]
	puts $out [string trim $new_time_out]
	close $out
}