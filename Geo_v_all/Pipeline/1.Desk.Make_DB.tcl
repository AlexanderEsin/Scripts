set direct /users/aesin/desktop

set org Geo_v_all/Input
set db_dir_forward Geo_fastas
set db_dir_reverse All_fastas

###########################################################################

cd $direct/$org
file mkdir DB_forward
file mkdir DB_reverse

cd $direct/$org/$db_dir_forward

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

cd $direct/$org/$db_dir_reverse

set dblist [glob *faa]
set x 1

foreach db $dblist {
	regsub -all " " $db "_" new_name
	file rename -force $db $new_name
	catch {exec makeblastdb -in $new_name -dbtype prot -out $direct/$org/DB_reverse/$db.db}
	puts "$new_name ==== $x"
	incr x
}

set x [expr $x -1]

puts "===DONE===\nDB_forward: $i sequences\nDB_reverse: $x sequences\n"