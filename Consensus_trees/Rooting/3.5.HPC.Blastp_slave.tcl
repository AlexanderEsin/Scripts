## Blastp all the proteomes against each geobacillus, perform the forward and reverse in the same script (requires only one splitting) ##
###########################################################################
set direct /users/aesin/desktop
set org Consensus_trees/Rooting/[lindex $argv 1]

set query_dir Fastas_forward
set query_dir_2 Fastas_reverse
set db_dir DB_split
set db_dir_2 DB_reverse
###########################################################################

cd $direct/$org
file mkdir Out_all
file mkdir Timing

set fold_no [lindex $argv 0]

cd $direct/$org/Input/$query_dir
set glist [glob *.faa]
puts $glist

cd $direct/$org/Input/$db_dir/$fold_no
set dblist [glob *.pin]
puts "$fold_no === [llength $dblist]"

foreach g $glist {
	foreach db $dblist {
		regsub {.db.pin} $db {} db
		
		if {$g != $db} {
			set timing_out ""

			regsub {.faa} $g {} x
			regsub {.faa} $db {} y
			set timer_forward [time {
				catch {exec blastp -query $direct/$org/Input/$query_dir/$g -db $direct/$org/Input/$db_dir/$fold_no/$db.db -out $direct/$org/Out_all/$x\&$y.tsv -evalue 1e-10 -outfmt 7 -max_target_seqs 1 -max_hsps 1 -num_threads 1}

				catch {exec blastp -query $direct/$org/Input/$query_dir_2/$db -db $direct/$org/Input/$db_dir_2/$g.db -out $direct/$org/Out_all/$y\&$x.tsv -evalue 1e-10 -outfmt 7 -max_target_seqs 1 -max_hsps 1 -num_threads 1}
			}]

			set time_forward [expr round([expr [lindex $timer_forward 0]/pow(10,6)])]

			set out [open $direct/$org/Timing/$db\.txt a]
			puts $out "$g\t$time_forward"
			close $out
		}
	}
}
