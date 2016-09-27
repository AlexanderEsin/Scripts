#!/bin/sh
#PBS -l walltime=70:00:00
#PBS -l select=1:ncpus=32:mem=80gb
#PBS -l place=free:shared:group=boardpair

##########################
direct=/home/ade110/Scripts
org=Geo_v_all
script=7.5.MCL_groups_slave.tcl
tclsh=1
bash=0
##########################

module load mcl/14.137

chmod +x $direct/$org/$script
# Make sure there is a space between each bit #
if [ "$tclsh" = 1 ]; then
	tclsh8.5 $direct/$org/$script
elif [ "$bash" = 1] ; then
	$direct/$org/$script
fi

