#!/bin/sh
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=1:mem=10gb

cd /mnt/storage/home/aesin/Work/Bacillus/Mowgli/Mowgli_output
# cd /mnt/storage/rawdata/Warnecke_Rawdata/ADE/Work/Staph

tar -cvzf Output_3.tar.gz Output_3