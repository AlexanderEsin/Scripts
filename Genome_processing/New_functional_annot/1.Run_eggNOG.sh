# Double check the ports for server and client

master="/Users/aesin/Desktop/Geo_analysis/Geo_ortholog_nucl/Functional_annot"
input_dir="Input_groups"
output_dir="Output_annots"

input_files=$(find $master/$input_dir -name "*faa")

for input in $input_files
do
	output_name=$(basename $input | sed 's/.faa/.txt/g')
	echo "Working on $output_name"

	python2.7 /users/aesin/Desktop/emapper/emapper.py -i $input \
	--output /Users/aesin/Desktop/Geo_analysis/Geo_ortholog_nucl/New_functional_annot/Output_annots/$output_name \
	-d /Users/aesin/Desktop/emapper/data/hmmdb_levels/bact_50/bact_50.hmm:localhost:53000 \
	--override \
	--cpu 18

	# wait 1000
done
