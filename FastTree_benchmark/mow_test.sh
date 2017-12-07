#!/bin/sh
#PBS -l walltime=1:00:00
#PBS -l select=1:ncpus=1:mem=1gb

module load mowgli/v2.0A

mowgli -h > $SCRATCH/mowgli_test.txt

