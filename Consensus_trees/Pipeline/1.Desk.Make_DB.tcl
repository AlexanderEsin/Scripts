#!/usr/local/bin/tclsh
set direct /users/aesin/desktop
set org Inside_group

set org Consensus_trees/$org/Input
set db_dir Fastas

###########################################################################

cd $direct/$org
file mkdir DB

cd $direct/$org/$db_dir

set dblist [glob *faa]
set i 1

foreach db $dblist {
	regsub -all " " $db "_" new_name
	file rename -force $db $new_name
	catch {exec makeblastdb -in $new_name -dbtype prot -out $direct/$org/DB/$new_name.db}
	puts "$new_name ==== $i"
	incr i
}

set i [expr $i -1]

puts "===DONE===\nDB: $i sequences\n"

###########################################################################

