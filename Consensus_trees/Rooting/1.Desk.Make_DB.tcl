#!/usr/local/bin/tclsh
set direct /users/aesin/desktop
set org [lindex $argv 0]

set org Consensus_trees/Rooting/$org/Input
set db_dir_forward Fastas_forward
set db_dir_reverse Fastas_reverse

###########################################################################

cd $direct/$org
file mkdir DB_forward
file mkdir DB_reverse

###########################################################################

cd $direct/$org/$db_dir_forward

set dblist [glob *faa]
set i 1

foreach db $dblist {
	regsub -all " " $db "_" new_name
	file rename -force $db $new_name
	catch {exec makeblastdb -in $new_name -dbtype prot -out $direct/$org/DB_reverse/$db.db}
	puts "$new_name ==== $i"
	incr i
}

set i [expr $i -1]

puts "===DONE===\nDB: $i sequences\n"

###########################################################################

cd $direct/$org/$db_dir_reverse

set dblist [glob *faa]
set i 1

foreach db $dblist {
	regsub -all " " $db "_" new_name
	file rename -force $db $new_name
	catch {exec makeblastdb -in $new_name -dbtype prot -out $direct/$org/DB_forward/$db.db}
	puts "$new_name ==== $i"
	incr i
}

set i [expr $i -1]

puts "===DONE===\nDB: $i sequences\n"

###########################################################################