#!/bin/sh

## I'm not entirely sure how many CPUs the DB builder can use (it ran at night, so I didn't get to see), but this seems like a safe bet re: cores and RAM 
ulimit -n 2048
##########################
direct=/users/aesin/desktop/Deer/Assembly

cd $direct
mkdir DB

for file in ./*fastq ; do
  if [ -e "$file" ] ; then   # Check whether file exists.

  	# Some file name manipulation ... #
  	file_name=$(echo $file | sed 's|./||')
  	#file_path="$(readlink -f $file_name)"
  	db_name=$(echo $file_name | sed 's|.fastq|_atram_db|')

  	# Make our own log file as the DB-making script doesn't have an inbuilt option like the main program #
  	echo "Now making database of $file_name ..." >> $direct/Log_db_make.txt
  	
  	/users/aesin/Desktop/Scripts/Deer/aTRAM-master/format_sra.pl -input $file_name -out DB/$db_name -max_processes 8 >> $direct/Log_db_make.txt 

  fi
done

echo "Database has been created on ..." >> $direct/Log_db_make.txt
date >> $direct/Log_db_make.txt