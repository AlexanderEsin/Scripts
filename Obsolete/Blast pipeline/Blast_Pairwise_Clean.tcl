cd ~/Desktop/test_blast/Input/marked_faa
set dblist [glob *.faa]
set i 1

set z [glob -nocomplain -type d *]
if {[regexp DB $z] == 1} {
} else {
	file mkdir DB
}
if {[regexp Out $z] == 1} {
} else {
	file mkdir Out
}

foreach db $dblist {
	exec makeblastdb -in $db -dbtype prot -out DB/$db.db
	puts "$i/[llength $dblist]"
	set glist [glob *.faa]
	foreach g $glist {
		if {$g == $db} {
			puts "skip"
		} else {
			regsub _protein.faa $g {} x
			regsub _protein.faa $db {} y
			puts "$x\_vs_$y"
			catch {exec blastp -query $g -db DB/$db.db -out Out/$x&$y.tsv -evalue 1e-10 -outfmt 7 -max_target_seqs 1 -max_hsps 1 -num_threads 3}
		}
	}
	incr i
}

file copy -force -- ~/Desktop/Test_Blast/Input/marked_faa/Out ~/Desktop/Test_Blast/Input
file delete -force -- ~/Desktop/Test_Blast/marked_faa/Out
file delete -force -- ~/Desktop/test_blast/Input/marked_faa/DB

puts "DONE"

