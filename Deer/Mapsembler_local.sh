#!/bin/sh

##########################
direct=/users/aesin/desktop/Deer/Assembly
split_no=$1
input_dir=Alpha_input_split/$split_no
reads_file="Ovirg_paired_interleaved.fastq"

cd $direct
output_dir=Mapsem_alpha_out/$split_no
mkdir -p $output_dir
mkdir -p Done_reads

cd $input_dir

for file in ./*.fasta ; do
  if [ -e "$file" ] ; then   # Check whether file exists.


  	# Some file name manipulation ... #
  	file_name=$(echo $file | sed 's|./||')
  	output_name=$(echo $file_name | sed 's|.fasta||')
  	echo "Now making assembly based on $file_path ..."

    cd $direct/$output_dir

  	/users/aesin/downloads/mapsembler2_pipeline/run_mapsembler2_pipeline.sh -s $direct/$input_dir/$file -r $direct/$reads_file -t 1 -p $output_name >> $direct/$output_dir/Log_assemble_short.txt

  	mv $direct/$input_dir/$file $direct/Done_reads/
  fi
done
