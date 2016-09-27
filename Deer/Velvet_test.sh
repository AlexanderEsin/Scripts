#!/bin/sh
#PBS -l walltime=00:10:00
#PBS -l select=1:ncpus=8:mem=10gb

###
module load velvet/1.2.10

echo $(module list)
echo $(velvetg)

exit
