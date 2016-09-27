set direct /users/aesin/desktop/Geobac/Outgroups/1_out_subtilis

proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}

## Bootstrap number ##

set N 1000

###########################################################################

cd $direct
file mkdir Tree_Build/$N
file mkdir Final_trees

cd $direct/Stitched_core_genome

set fams [glob *.fas]

foreach fam $fams {
	regsub ".fas" $fam {.txt} fname
	catch {exec raxml -f a -s $fam -n $fname -m PROTCATAUTO -p 1234 -x 10001 -N $N -T 7}
}

set reduced [glob -nocomplain *.reduced]
foreach red $reduced {
	file delete $red
}

set tree_files [glob RAxML*]
foreach file $tree_files {
	file rename $file $direct/Tree_Build/$N/$file
}

set key_file [lindex [glob *KEY*] 0]
openfile $key_file
set key [string trim $data]
file copy -force $key_file $direct/Tree_Build/$N/

cd $direct/Tree_Build/$N

set pata {k+?[0-9]+?:+?}
set patb {\t+.+}

set besttrees [glob -nocomplain *bipartitions.*]
if {[llength $besttrees] > 0} {
	foreach tree $besttrees {
		openfile $tree
		set branchids [regexp -all -inline $pata $data]
		foreach k_val $branchids {
			regsub {:} $k_val {} k_val_trunc
			regexp -line $k_val_trunc$patb $key hit
			set binomial [lindex [split [string trim $hit] \t] 1]
			set id [lindex [split [string trim $hit] \t] 2]
			# If you use k_val_trunc below, you may substitute k17 when you want to find k1 - so searching for k1: is more accurate #
			regsub $k_val $data "$binomial\:" data
		}
		set out [open Geobacillus_core_tree_N$N.txt w]
		puts $out $data
		close $out
		file copy -force Geobacillus_core_tree_N$N.txt $direct/Final_trees/
	}
}
