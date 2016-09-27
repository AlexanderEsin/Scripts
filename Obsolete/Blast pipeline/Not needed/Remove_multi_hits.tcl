proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}

cd ~/
cd ~/desktop/test_blast/input/out/pairwise_7_names
set filelist [glob *.tsv]

foreach tsv $filelist {
	openfile $tsv
	set genelist ""
	set out ""
	set lineparse [split $data "\n"]
	regsub "{}" $lineparse "" lineparse
	foreach line $lineparse {
		set x [lindex $line 0]
		if {$x == "#"} {
			append out "$line\n"
		} elseif {[regexp $x $genelist] == 1} {
			regsub -all $line {} na
		} else {
			set gene [lindex $line 0]
			append genelist "$gene\t"
			append out "$line\n"
		}
	}
	set out2 [open ~/Desktop/Test_Blast/Input/Out/Pairwise_7_names_1hit/$tsv a]
	puts $out2 "$out"
	close $out2
}