## Blastp all the proteomes against each geobacillus, perform the forward and reverse in the same script (requires only one splitting) ##
###########################################################################
set direct /users/aesin/desktop/Geo_v_all/2.0
set org Geo_only_groups; #Anoxy_geo_only_groups or Geo_only_groups
set blast_v_what Blast_v_Gthermo

set org $org/$blast_v_what

set query_forward_dir Split
set db_forward_dir DB

set fold_no [lindex $argv 0]

###########################################################################

cd $direct/$org

puts "$fold_no == [exec pwd]"

file mkdir Out_all
file mkdir Timing

cd $direct/$org/Input/$query_forward_dir/$fold_no
set glist [glob *.faa]

cd $direct/$org/Input/$db_forward_dir
set dblist [glob *.pin]

foreach g $glist {
	foreach db $dblist {
		regsub {.pin} $db {} db

		if {$g != $db} {
			set timing_out ""

			regsub {.faa} $g {} x
			regsub {.faa} $db {} y
			set timer_forward [time {
				catch {exec blastp -query $direct/$org/Input/$query_forward_dir/$fold_no/$g -db $direct/$org/Input/$db_forward_dir/$db -out $direct/$org/Out_all/$x\&$y.tsv -evalue 1e-05 -outfmt 7 -max_target_seqs 1 -max_hsps 1 -num_threads 1}
			}]

			set time_forward [expr round([expr [lindex $timer_forward 0]/pow(10,6)])]

			set out [open $direct/$org/Timing/$db\.txt a]
			puts $out "$g\t$time_forward"
			close $out
		}
	}
}

puts "\nFolder: $fold_no == DONE"