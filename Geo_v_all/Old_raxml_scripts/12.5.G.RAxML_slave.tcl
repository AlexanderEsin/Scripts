set eval -100
set ival 2.5
set direct /scratch/ade110/Geo_v_all/$eval
set size giant; #small, large or giant

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

set start_time [clock seconds]
set soft_end_time [expr round($start_time + (71*60*60))]

cd $direct
file mkdir Trees_$ival/$size
file mkdir Tree_build_$ival/$size

cd $direct/MSA_split_$ival/$size/[lindex $argv 0]

set dirs [glob -nocomplain -type d *]
if {[llength $dirs] == 0} {
	exit
}

set dirs [lsort -dictionary $dirs]
reverse_dict $dirs
set dirs $revdict

set current_run_start 0
set current_run_finish 0

foreach directory $dirs {
	set current_time [clock seconds]
	set previous_run_time [expr $current_run_finish - $current_run_start]
	set projected_run_time [expr round($previous_run_time * 3)]
	set projected_finish_time [expr $current_time + $projected_run_time]

	if {$projected_finish_time > $soft_end_time} {
		puts "DONE: $projected_finish_time > $soft_end_time ==== $current_time"
		exit 2
	}

	###

	set current_run_start [clock seconds]
	cd $directory
	set align [lindex [glob *.fas] 0]
	regsub {.fas} $align {} out_number
	catch {exec raxml -f d -s $align -n $out_number\.txt -m PROTCATAUTO -p 1234 -T 16}

	set key [lindex [glob KEY*] 0]
	openfile $key
	set key $data

	set pata {k+?[0-9]+?:+?}
	set patb {\t+?.+?\n+?}

	set besttrees [glob -nocomplain *bestTree.*]
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
				regsub $k_val $data "$binomial \{$id\}:" data
			}
			set out [open $out_number\_tree.txt w]
			puts $out $data
			close $out
			file copy $out_number\_tree.txt $direct/Trees_$ival/$size/
		}
	}	
	cd ..
	file rename -force $directory $direct/Tree_build_$ival/$size/
	set current_run_finish [clock seconds]
}


