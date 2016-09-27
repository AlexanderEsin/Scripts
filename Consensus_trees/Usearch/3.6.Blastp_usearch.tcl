## Blastp all the proteomes against each geobacillus, perform the forward and reverse in the same script (requires only one splitting) ##
###########################################################################
set direct /users/aesin/Desktop
set org Consensus_trees/Bacilli_u

set query_dir Split
set db_dir DB
###########################################################################

cd $direct/$org
file mkdir Out_all
file mkdir Timing

set fold_no [lindex $argv 0]

cd $direct/$org/Input/$query_dir/$fold_no
set glist [glob *.faa]

cd $direct/$org/Input/$db_dir
set dblist [glob *.udb]

foreach g $glist {
	foreach db $dblist {
		regsub {.udb} $db {} db
		
		if {$g != $db} {
			set timing_out ""

			regsub {.faa} $g {} x
			regsub {.faa} $db {} y
			set timer_forward [time {
				catch {exec usearch -usearch_local $direct/$org/Input/$query_dir/$fold_no/$g -db $direct/$org/Input/$db_dir/$db.udb -id 0.25 -evalue 1e-10 -blast6out $direct/$org/Out_all/$x\&$y.tsv}
			}]

			set time_forward [expr round([expr [lindex $timer_forward 0]/pow(10,6)])]

			set out [open $direct/$org/Timing/$db\.txt a]
			puts $out "$g\t$time_forward"
			close $out
		}
	}
}