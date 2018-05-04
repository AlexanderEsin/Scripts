## Run the progressiveMauve command line. The binary is run out of the app directory - added to $PATH ##

cd "/Users/aesin/Desktop/Geo_again/Genomes/AG_genome_noPlasmid_gbffs"
genome_files=*

for n in $genome_files; do
	new_name=$(echo $n | sed 's/.gbff/.gbk/')
	mv $n $new_name
done

ulimit -n 10000

## Remember to remove any .sslist files on reruns
progressiveMauve --output=../AG_mauveAlign/AG_25_mauveNoPlasmid --disable-cache $genome_files

