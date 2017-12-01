#!/bin/bash
# Create blast DB for each of the AnoGeo proteomes
direct=/users/aesin/Desktop/Geo_again/BLAST
mkdir -p $direct/All_DB
in_files="$(find $direct/IN_all -name "*.fasta" -print0 | xargs -0 ls)"
for input_fasta in $in_files
do
	db_name=$(basename $(echo $input_fasta | sed 's/_proteome.fasta//g'))

	if [ -e $direct/All_DB/$db_name\.phr ]; then
		echo "DB exists..."
		continue
	else
		makeblastdb -in $input_fasta -dbtype prot -out $direct/All_DB/$db_name 
	fi
	
done