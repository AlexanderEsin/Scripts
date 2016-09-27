set direct /users/aesin/desktop/Geobac/Geobac_1out_difference

proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}

###########################################################################
#Blastp the genes that are no longer picked out as orthologous when the first outgroup/
#is included in the initial screen. Then screen for any significant blast hits to the/
#1st outgroup and pull out those that don't have apparent homology into "No_homologs"
###########################################################################

cd $direct

set z [glob -nocomplain -type d *]
if {[regexp DB $z] == 1} {
} else {
	file mkdir DB
}
if {[regexp Out $z] == 1} {
} else {
	file mkdir Out
}

cd $direct/DB_fasta
set dblist [glob *.faa]
set i 1

foreach db $dblist {
	exec makeblastdb -in $db -dbtype prot -out $direct/DB/$db.db
	puts "$i/[llength $dblist]"
	cd $direct/Group_fastas
	set glist [glob *.txt]
	foreach g $glist {
		if {$g == $db} {
			puts "skip"
		} else {
			regsub {.txt} $g {} x
			set y "Bacillus_smithii"
			puts "$x\_vs_$y"
			catch {exec blastp -query $g -db $direct/DB/$db.db -out $direct/Out/$x&$y.tsv -evalue 1e-10 -outfmt 7 -max_target_seqs 1 -max_hsps 1 -num_threads 7}
		}
	}
	incr i
}

file delete -force -- $direct/DB

puts "DONE"

###########################################################################

cd $direct

file mkdir No_homologs

cd $direct/Out

set blast_results [glob *.tsv]

foreach result $blast_results {
	openfile $result
	set x [regexp -all "0 hits found" $data]
	if {$x == 16} {
		regsub {&Bacillus_smithii.tsv} $result {.txt} orig_fasta
		file copy $direct/Group_fastas/$orig_fasta $direct/No_homologs/
	}
}