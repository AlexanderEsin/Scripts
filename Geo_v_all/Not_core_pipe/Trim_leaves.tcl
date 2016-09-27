proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}

set pata {\s+?\{+?.+?\}+?}

cd ~/desktop/Prun_test
openfile [lindex [glob *tree*] 0]

set data [string trim $data]

set test [regexp -all -inline $pata $data]

regsub -all $pata $data {} leaf_trim
puts $leaf_trim

set out [open ~/desktop/Prun_test/Leaf_trim_test_tree.txt w]
puts $out $leaf_trim
close $out