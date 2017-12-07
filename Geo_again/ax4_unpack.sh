#!/bin/sh
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=10gb

cd /scratch2/ade110/[OLD]/Geo_v_all
tar -cvzf Out_all.tar.gz Out_all/
