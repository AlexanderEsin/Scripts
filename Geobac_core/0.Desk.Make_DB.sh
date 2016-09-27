#!/bin/sh

org=Geobac/All_geobac/Input
db_dir_forward=Fastas
#db_dir_reverse=All_fastas

###########################################################################
### Make both the forward and (if needed) reverse blast databases ###

cd ~/desktop/$org
mkdir -p DB_forward
#mkdir -p DB_reverse

cd ~/desktop/$org/$db_dir_forward

dblist=?*.faa
i=1

for db in $dblist
do
	# Remove all spaces from the file name #
	new_name=$(echo $db | sed 's| |_|g')
	mv "$db" "$new_name"
	makeblastdb -in $new_name -dbtype prot -out ~/desktop/$org/DB_forward/$db.db
	printf "$new_name ==== $i"
	i=$((i+1))
done
i=$((i-1))

#cd ~/desktop/$org/$db_dir_reverse

#dblist=?*.faa
#x=1

#for db in $dblist
#do
	# Remove all spaces from the file name #
#	new_name=$(echo $db | sed 's| |_|g')
#	mv "$db" "$new_name"
#	makeblastdb -in $new_name -dbtype prot -out ~/desktop/$org/DB_reverse/$db.db
#	printf "$new_name ==== $x"
#	x=$((x+1))
#done
#x=$((x-1))


#printf "===DONE===\nDB_forward: $i sequences\nDB_reverse: $x sequences\n"


