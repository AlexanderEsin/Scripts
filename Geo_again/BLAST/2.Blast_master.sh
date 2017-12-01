#!/bin/sh
#PBS -l walltime=3:00:00
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -J 1-1691

module load blast+/2.2.28

chmod +x /home/ade110/Scripts/Geo_again/BLAST/2.Blast_slave.tcl

new_number=$(printf %04d $PBS_ARRAY_INDEX)
echo $new_number

/home/ade110/Scripts/Geo_again/BLAST/2.Blast_slave.tcl $new_number