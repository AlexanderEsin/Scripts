#!/bin/sh
#PBS -l walltime=1:00:00
#PBS -l select=1:ncpus=1:mem=5gb


cd /scratch/ade110/Geo_again/BLAST
tar -xvzf DB_all.tar.gz
