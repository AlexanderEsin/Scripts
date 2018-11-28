#!/bin/bash
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=48:mem=40gb

module load raxml

master_dir="/home/ade110/Work/alignTW"
#master_dir="/Users/aesin/Desktop/alignTW"
input_dir="input"
output_dir="output"

# Make output directory
mkdir -p $master_dir/$output_dir

# List of alignment files
align_files=$(find $master_dir/$input_dir -name "*doublets*.afa" -print0 | xargs -0 ls)
echo $align_files


for align in $align_files
do
	fileName=$(basename $align)
	noExt="${fileName%.*}"

	raxmlHPC-PTHREADS-AVX2 -f a -s $align -n $noExt -w $master_dir/$output_dir -m PROTGAMMAAUTO -p $RANDOM -x $RANDOM -N 100 -T 20
done

echo "DONE"
