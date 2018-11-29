#!/bin/sh
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=2:mem=20gb
#PBS -J 1-2

#6724
module load fasttree

chmod +x /home/ade110/Scripts/Geo_again/Align_tree/9.Group_fastTree_slave.tcl

/home/ade110/Scripts/Geo_again/Align_tree/9.Group_fastTree_slave.tcl $PBS_ARRAY_INDEX