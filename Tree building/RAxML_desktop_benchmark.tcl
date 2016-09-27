set direct /users/aesin/desktop/

set size small; #large or small

proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}

proc reverse_dict {lst} {
	global revdict
	while {[llength $lst] > 0} {
		set element [lindex $lst end]
		lappend revdict $element
		set lst [lrange $lst 0 end-1]
	}
}

###########################################################################

cd $direct
file mkdir Final_trees/small/
file mkdir Tree_build

cd $direct/MSA_split/$size/1

set dirs [glob -type d *]
set dirs [lsort -dictionary $dirs]
reverse_dict $dirs
set dirs $revdict

set dirs {1085}

foreach directory $dirs {
	cd $directory
	set align [lindex [glob *.fas] 0]
	regsub {.fas} $align {} out_number
	catch {exec raxml -f a -s $align -n $out_number\.txt -m PROTCATDAYHOFF -p 1234 -x 10001 -N 100 -T 4}

	set key [lindex [glob KEY*] 0]
	openfile $key
	set key $data

	set pata {k+?[0-9]+?:+?}
	set patb {\t+?.+?\n+?}

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
				regsub $k_val_trunc $data "$binomial \{$id\}" data
			}
			set out [open $out_number\_tree.txt w]
			puts $out $data
			close $out
			file copy $out_number\_tree.txt $direct/Final_trees/
		}
	}	
	cd ..
	file rename -force $directory $direct/Tree_build
}

