set ival 2.0
set direct /scratch/ade110/Geo_v_all/$ival/MSA_split


proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}
###########################################################################

cd $direct
set core_splits [glob -nocomplain -type d *core*]
puts "REMAINING ALIGNMENTS TO BE RECONSTRUCTED:"

foreach core_split $core_splits {
	cd $direct/$core_split
	set sizes [glob -type d *]

	foreach size $sizes {
		cd $direct/$core_split/$size

		set dirs [glob -nocomplain -type d *]
		set total_no 0
		foreach dir $dirs {
			cd $dir
			set file_no [llength [glob -nocomplain -type d *]]
			if {$file_no != 0 && $size eq "small"} {
				puts "SMALL: $dir == $file_no"
			}
			set total_no [expr $file_no + $total_no]
			cd ..
		}

		puts "$core_split === $size === $total_no"
	}
	puts "\n"
}