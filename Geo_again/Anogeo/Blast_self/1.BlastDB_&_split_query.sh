#!/bin/bash

# Functions
ceiling_divide () {
  ceiling_result=$((($1+$2-1)/$2))
  echo $ceiling_result
}


# Create blast DB for each of the AnoGeo proteomes
direct=/users/aesin/Desktop/Geo_again/Anogeo_analysis/Blast_anogeo
mkdir -p $direct/DB

in_files="$(find $direct/Input -name "*.fasta" -print0 | xargs -0 ls)"
for input_fasta in $in_files
do
	db_name=$(basename $(echo $input_fasta | sed 's/_proteome.fasta//g'))
	makeblastdb -in $input_fasta -dbtype prot -out $direct/DB/$db_name
done


# Split input queries into multiple subfolders for crude parallelisation
num_subdir=20
mkdir -p $direct/Input_split

file_total=$(echo $in_files | wc -w)
per_folder=$(ceiling_divide $file_total $num_subdir)

# Split the files into multiple sub-directories
i=0
for f in $in_files
do
    d=dir_$(printf %03d $((i/$per_folder+1)))
    mkdir -p $direct/Input_split/$d
    cp "$f" $direct/Input_split/$d
    let i++
done


