#!/bin/bash
# genbank_to_fasta.py only works with biopython version 1.67 - installed to hombrew python2
direct=/users/aesin/Desktop/Geo_again/Genomes/Anogeo
in_files="$(find $direct/Anogeo_genomes_raw/Genomic_gbff -name "*.gbff.gz" -print0 | xargs -0 ls)"

for in_file in $in_files
do
	gunzip $in_file
	in_file=$(echo $in_file | sed 's/_genomic.gbff.gz/_genomic.gbff/g')
	out_file=$(echo $in_file | sed 's/_genomic.gbff/_proteome.fasta/g')
	genbank_to_fasta.py -i $in_file -o $out_file -m 'genbank' -s 'aa' -f 'CDS' -d 'spacepipe' -q 'protein_id,locus_tag,gene,product,location'
done

mv $direct/Anogeo_genomes_raw/Genomic_gbff/*fasta $direct/Anogeo_proteomes/Proteome_fasta_raw