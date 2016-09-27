## Blastp all the proteomes against each geobacillus, perform the forward and reverse in the same script (requires only one splitting) ##
###########################################################################
set direct /scratch/ade110
set org Geo_v_all

set query_forward_dir Split
set db_forward_dir DB_forward

set query_reverse_dir Geo_fastas
set db_reverse_dir DB_reverse

###########################################################################

cd $direct/$org
file mkdir Out_all
file mkdir Timing

set fold_no [lindex $argv 0]

cd $direct/$org/Input/$query_forward_dir/$fold_no
set glist [glob *.faa]

cd $direct/$org/Input/$db_forward_dir
set dblist [glob *.pin]


foreach g $glist {
	foreach db $dblist {
		regsub {.db.pin} $db {} db
		
		if {$g != $db} {
			set timing_out ""

			regsub {.faa} $g {} x
			regsub {.faa} $db {} y
			set timer_forward [time {
				catch {exec blastp -query $direct/$org/Input/$query_forward_dir/$fold_no/$g -db $direct/$org/Input/$db_forward_dir/$db.db -out $direct/$org/Out_all/$x\&$y.tsv -evalue 1e-10 -outfmt 7 -max_target_seqs 1  -max_hsps_per_subject 1 -seg yes -soft_masking true -num_threads 1}
			}]

			set time_forward [expr round([expr [lindex $timer_forward 0]/pow(10,6)])]
		

			set timer_reverse [time {
				catch {exec blastp -query $direct/$org/Input/$query_reverse_dir/$db -db $direct/$org/Input/$db_reverse_dir/$g.db -out $direct/$org/Out_all/$y\&$x.tsv -evalue 1e-10 -outfmt 7 -max_target_seqs 1 -max_hsps_per_subject 1 -seg yes -soft_masking true -num_threads 1}
			}]

			set time_reverse [expr round([expr [lindex $timer_reverse 0]/pow(10,6)])]

			set out [open $direct/$org/Timing/$db\.txt a]
			puts $out "$g\t$time_forward\t$time_reverse"
			close $out
		}
	}
}