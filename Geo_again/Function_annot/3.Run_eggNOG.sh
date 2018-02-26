# Double check the ports for server and client

master="/Users/aesin/Desktop/Geo_again/Functional_annotation"
input_dir="AG_only_groups"
output_dir="Output_annots"
mkdir -p $master/$output_dir

input_files=$(find $master/$input_dir -name "*.fasta")

for input in $input_files
do
	output_name=$(basename $input | sed 's/.fasta/.txt/g')
	echo "Working on $output_name"

	python2.7 /users/aesin/Desktop/emapper/emapper.py -i $input \
	--output $master/$output_dir/$output_name \
	-d /Users/aesin/Desktop/emapper/data/hmmdb_levels/bact_50/bact_50.hmm:localhost:51500 \
	--override

	# wait 1000
done
