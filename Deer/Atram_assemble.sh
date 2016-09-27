#!/bin/sh
#PBS -l walltime=70:00:00
#PBS -l select=1:ncpus=16:ompthreads=16:mem=800gb

module load atram/1.0.4
module switch velvet/1.1.04 velvet/1.2.10

##########################
direct=/scratch/ade110/Deer

cd $direct/DB
## Get the database variable ##
db_name=$(ls | sort -n | head -1 | sed 's|\..*||')

cd $direct
output_dir=Output_alpha_OMP
mkdir $output_dir

for file in ./*alpha_1.fasta ; do
  if [ -e "$file" ] ; then   # Check whether file exists.

  	# Some file name manipulation ... #
  	file_name=$(echo $file | sed 's|./||')
  	output_name=$(echo $file_name | sed 's|.fasta||')
  	file_path="$(readlink -f $file_name)"
  	echo "Now making assembly based on $file_path ..."

  	## Set max_memory (in GBs) to 20GB below requested max - some users reported spiking above max argument given. Default assembler is Velvet (documentation), so no need to specify ##
  	aTRAM_AE.pl -reads $direct/DB/$db_name -target $file_name -kmer 21 -ins_length 270 -max_target_seqs 10000 -exp_coverage 9 -max_processes 16 -max_memory 780 -iterations 4 -output $direct/$output_dir/Out_atram_$output_name -log_file $direct/$output_dir/Log_$file_name\_assembly.txt  >> $direct/$output_dir/Log_assemble_short.txt
  fi
done
