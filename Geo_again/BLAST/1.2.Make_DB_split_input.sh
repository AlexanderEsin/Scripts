#!/bin/bash

## Use the short queue to run the Geo_again blast

# Functions
ceiling_divide () {
  ceiling_result=$((($1+$2-1)/$2))
  echo $ceiling_result
}


# Create blast DB for each of the AnoGeo proteomes
direct=/users/aesin/Desktop/Geo_again/BLAST
mkdir -p $direct/DB_all
in_files="$(find $direct/IN_all -name "*.fasta" -print0 | xargs -0 ls)"

## For each fasta check if DB already exists
## If not, create blast DB
for input_fasta in $in_files
do
	db_name=$(basename $(echo $input_fasta | sed 's/_proteome.fasta//g'))

	if [ -e $direct/DB_all/$db_name\.phr ]; then
		echo "DB exists..."
		continue
	else
		makeblastdb -in $input_fasta -dbtype prot -out $direct/DB_all/$db_name 
	fi
	
done



# Split input queries into multiple subfolders for crude parallelisation
total_in=$(echo $in_files | wc -w)
num_subdir=$(ceiling_divide $total_in 3)

per_folder=3

# Split the files into multiple sub-directories
i=0
for f in $in_files
do
    d=dir_$(printf %04d $((i/$per_folder+1)))
    mkdir -p $direct/IN_all_split/$d
    cp "$f" $direct/IN_all_split/$d
    let i++
done


