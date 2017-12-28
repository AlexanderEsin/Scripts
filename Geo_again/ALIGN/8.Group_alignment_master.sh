#!/bin/sh
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=4:mem=20gb
#PBS -J 1-6726

module load clustal-omega/1.2

chmod +x /home/ade110/Scripts/Geo_again/ALIGN/8.Group_alignment_slave.tcl

/home/ade110/Scripts/Geo_again/ALIGN/8.Group_alignment_slave.tcl $PBS_ARRAY_INDEX