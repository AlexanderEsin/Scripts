## Blastp all the proteomes against each other ##
###########################################################################
set direct /csc/rawdata/Warnecke/ADE
set org Consensus_trees/Inside_group

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
set dblist [glob *.pin]

foreach g $glist {
	foreach db $dblist {
		regsub {.db.pin} $db {} db
		
		if {$g != $db} {
			set timing_out ""

			regsub {.faa} $g {} x
			regsub {.faa} $db {} y
			set timer_forward [time {
				catch {exec blastp -query $direct/$org/Input/$query_dir/$fold_no/$g -db $direct/$org/Input/$db_dir/$db.db -out $direct/$org/Out_all/$x\&$y.tsv -evalue 1e-10 -outfmt 7 -max_target_seqs 1 -max_hsps_per_subject 1 -num_threads 1}
			}]

			set time_forward [expr round([expr [lindex $timer_forward 0]/pow(10,6)])]

			set out [open $direct/$org/Timing/$db\.txt a]
			puts $out "$g\t$time_forward"
			close $out
		}
	}
}