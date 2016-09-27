set direct /users/aesin/desktop
set org Bacilli_u

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
	catch {exec usearch -makeudb_usearch $new_name -output $direct/$org/DB/$db.udb}
	puts "$new_name ==== $i"
	incr i
}

set i [expr $i -1]

puts "===DONE===\nDB: $i sequences\n"


# usearch -makeudb_usearch db.fasta -output db.udb