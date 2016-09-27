#!/bin/sh
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=40:mem=420gb

## Run trinity ##

module load bowtie/1.0.0
module load samtools/1.2
module load intel-suite/2015
module load java

cd /csc/analysis/Warnecke/PDD_rnaseq

mkdir -p Trinity_out

echo $(which samtools) > Trinity_out/run.log

/csc/analysis/Warnecke/trinityrnaseq/Trinity --seqType fq --left Paired_val/1_S1_L001_R1_001_val_1.fq --right Paired_val/1_S1_L001_R2_001_val_2.fq --output Trinity_out/ --CPU 40 --max_memory 420G >> Trinity_out/run.log 2>&1
