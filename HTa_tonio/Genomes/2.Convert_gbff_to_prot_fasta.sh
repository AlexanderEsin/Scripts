#!/bin/bash

## Set the input + output directories
org="protozoa"

direct="/users/aesin/Desktop/HTa_tonio"
genome_in_dir="$direct/Genomes/$org/Genome_gbffs"
proteome_out_dir="$direct/Proteomes/Proteome_fastas"

## Make the output directory, if it doesn't exist
mkdir -p $proteome_out_dir

## The list of input gbff files
in_gbff_files="$(find $genome_in_dir -name "*.gbff.gz" -print0 | xargs -0 ls)"

## For each file, unzip and convert to fasta with genbank_to_fasta.py
i=1
for in_file in $in_gbff_files; do
	
	## Prep input and output names
	unzip_name=$(echo $in_file | sed 's/_genomic.gbff.gz/_genomic.gbff/g')
	out_file_name=$(echo $unzip_name | sed 's/_genomic.gbff/_proteome.fasta/g')
	out_file_path="$proteome_out_dir/$(basename $out_file_name)"


	# if [ -e $out_file_path ]; then
	# 	echo "DONE"
	# 	continue
	# else
	echo "Converting GBFF to FASTA.... DONE: $i"
	## Unzip file
	gunzip $in_file

	## Convert to fasta
	genbank_to_fasta.py -i $unzip_name -o $out_file_path -m 'genbank' -s 'aa' -f 'CDS' -d 'spacepipe' -q 'protein_id,locus_tag,gene,product,location' &>/dev/null

	## Gzip the file back
	pigz $unzip_name
	let i++
	# fi
done

echo -e "\nALL DONE :)"



	