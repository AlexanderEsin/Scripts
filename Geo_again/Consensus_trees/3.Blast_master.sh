#!/bin/sh
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=4:mem=20gb
#PBS -J 1-43

module load blast+/2.2.28

group_name=$GROUP

chmod +x /home/ade110/Scripts/Geo_again/Consensus_trees/3.Blast_slave.tcl

new_number=$(printf %03d $PBS_ARRAY_INDEX)
echo $new_number

/home/ade110/Scripts/Geo_again/Consensus_trees/3.Blast_slave.tcl $group_name $new_number