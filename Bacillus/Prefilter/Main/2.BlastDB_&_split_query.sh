#!/bin/bash

# Functions
ceiling_divide () {
  ceiling_result=$((($1+$2-1)/$2))
  echo $ceiling_result
}


## Directories
direct=/users/aesin/Desktop/Bacillus/Prefilter/Bacillus_withOut
blast_in_dir=$direct/BLAST/IN_all
mkdir -p $blast_in_dir

## Copy all the clean proteomes to the BLAST directory
clean_fastas="$(find $direct/Proteomes/Clean -name "*.fasta" -print0 | xargs -0 ls)"
for fasta in $clean_fastas
do
	cp -f $fasta $blast_in_dir
done



## Create blast DB for each of group proteomes
in_files="$(find $blast_in_dir -name "*.fasta" -print0 | xargs -0 ls)"
blast_db_dir=$direct/BLAST/DB_all
mkdir -p $blast_db_dir

## For each fasta check if DB already exists
## If not, create blast DB
for input_fasta in $in_files
do
	db_name=$(basename $(echo $input_fasta | sed 's/_proteome.fasta//g'))

	if [ -e $blast_db_dir/$db_name\.phr ]; then
		echo "DB exists..."
		continue
	else
		makeblastdb -in $input_fasta -dbtype prot -out $blast_db_dir/$db_name 
	fi
done


# Split input queries into multiple subfolders for crude parallelisation
num_subdir=50
in_split_dir=$direct/BLAST/IN_all_split
mkdir -p $in_split_dir

file_total=$(echo $in_files | wc -w)
per_folder=$(ceiling_divide $file_total $num_subdir)

# Split the files into multiple sub-directories
i=0
for f in $in_files
do
    d=dir_$(printf %03d $((i/$per_folder+1)))
    mkdir -p $in_split_dir/$d
    cp "$f" $in_split_dir/$d
    let i++
done


