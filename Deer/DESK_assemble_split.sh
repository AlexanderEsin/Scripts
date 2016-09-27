##########################
direct=/users/aesin/desktop/Deer/Assembly
input_dir=Alpha_input_round2
db_dir=$direct/DB

cd $db_dir
## Get the database variable ##
db_name=$(ls | sort -n | head -1 | sed 's|\..*||')

cd $direct
output_dir=Alpha_output_round2
done_dir=Alpha_done_reads_round2
mkdir -p $output_dir
mkdir -p $done_dir

cd $input_dir

for file in ./*.fasta ; do
  if [ -e "$file" ] ; then   # Check whether file exists.

    # Some file name manipulation ... #
    file_name=$(echo $file | sed 's|./||')
    output_name=$(echo $file_name | sed 's|.fasta||')
    echo "Now making assembly based on $file_name ..."

    ## Set max_memory (in GBs) below requested max - some users reported spiking above max argument given. 
    ## Default assembler is Velvet (documentation), so no need to specify
    ## Kmer size == 21
    ## Insert length == 270
    /users/aesin/downloads/aTRAM-master/aTRAM_AE.pl -reads $db_dir/$db_name -target $file_name -kmer 31 -max_target_seqs 2000 -ins_length 270 -exp_coverage 8 -cov_cutoff 2 -output $direct/$output_dir/Out_atram_$output_name -log_file $direct/$output_dir/Log_$file_name\_assembly.txt -velvet_log $direct/$output_dir/Velvet_log_$file_name\.txt -max_memory 24 -max_processes 7 -iterations 5 >> $direct/$output_dir/Log_assemble_short.txt

    mv $file $direct/$done_dir
  fi
done

