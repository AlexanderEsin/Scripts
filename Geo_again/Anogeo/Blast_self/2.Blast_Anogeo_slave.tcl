#!/usr/local/bin/tclsh

## Blastp all the Geobacilli against one another ##
###########################################################################
set direct /users/aesin/Desktop/Geo_again/Anogeo_analysis/Blast_anogeo

set query_dir Input_split
set db_dir DB

###########################################################################

cd $direct
file mkdir Blast_out
file mkdir Timing

set fold_no [lindex $argv 0]

cd $direct/$query_dir/dir_$fold_no
set proteome_l [glob *.fasta]

cd $direct/$db_dir
set db_l [glob *.pin]


foreach query $proteome_l {
	set query_base [string range $query 0 end-15]

	foreach db $db_l {
		set db_base [string range $db 0 end-4]
		
		if {$query_base != $db_base} {

			set timing_out ""

			set timer_forward [time {
				catch {exec blastp -query $direct/$query_dir/dir_$fold_no/$query -db $direct/$db_dir/$db_base -out $direct/Blast_out/$query_base\&$db_base.tsv -evalue 1e-10 -outfmt 6 -max_target_seqs 1 -max_hsps 1 -num_threads 1}
			}]

			set time_forward [expr round([expr [lindex $timer_forward 0]/pow(10,6)])]
		
			set out [open $direct/Timing/$db_base\.txt a]
			puts $out "$query\t$time_forward"
			close $out
		}
	}
}