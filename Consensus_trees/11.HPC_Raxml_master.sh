#!/bin/sh
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=16:mem=16gb
#PBS -l place=free:shared:group=board

module load raxml/2015-02-17 

org=Consensus_trees

chmod +x /home/ade110/Scripts/$org/11.5.HPC_Raxml_slave.tcl

tclsh8.5 /home/ade110/Scripts/$org/11.5.HPC_Raxml_slave.tcl