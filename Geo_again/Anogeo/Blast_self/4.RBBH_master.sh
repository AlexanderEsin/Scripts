#!/bin/bash

chmod +x ~/Documents/Scripts/Geo_again/Anogeo/Blast_self/4.RBBH_calculate_slave.tcl

ival=$1
eval=$2
trunc_eval=$(echo ${eval:2})

echo "Calculating RBBH hits ..."

~/Documents/Scripts/Geo_again/Anogeo/Blast_self/4.RBBH_calculate_slave.tcl $eval 1 &
~/Documents/Scripts/Geo_again/Anogeo/Blast_self/4.RBBH_calculate_slave.tcl $eval 2 &
~/Documents/Scripts/Geo_again/Anogeo/Blast_self/4.RBBH_calculate_slave.tcl $eval 3 &
~/Documents/Scripts/Geo_again/Anogeo/Blast_self/4.RBBH_calculate_slave.tcl $eval 4 &
~/Documents/Scripts/Geo_again/Anogeo/Blast_self/4.RBBH_calculate_slave.tcl $eval 5 &
~/Documents/Scripts/Geo_again/Anogeo/Blast_self/4.RBBH_calculate_slave.tcl $eval 6 &
~/Documents/Scripts/Geo_again/Anogeo/Blast_self/4.RBBH_calculate_slave.tcl $eval 7 &
~/Documents/Scripts/Geo_again/Anogeo/Blast_self/4.RBBH_calculate_slave.tcl $eval 8 &
~/Documents/Scripts/Geo_again/Anogeo/Blast_self/4.RBBH_calculate_slave.tcl $eval 9 &
~/Documents/Scripts/Geo_again/Anogeo/Blast_self/4.RBBH_calculate_slave.tcl $eval 10 &
~/Documents/Scripts/Geo_again/Anogeo/Blast_self/4.RBBH_calculate_slave.tcl $eval 11 &
~/Documents/Scripts/Geo_again/Anogeo/Blast_self/4.RBBH_calculate_slave.tcl $eval 12 &
~/Documents/Scripts/Geo_again/Anogeo/Blast_self/4.RBBH_calculate_slave.tcl $eval 13 &
~/Documents/Scripts/Geo_again/Anogeo/Blast_self/4.RBBH_calculate_slave.tcl $eval 14 &
~/Documents/Scripts/Geo_again/Anogeo/Blast_self/4.RBBH_calculate_slave.tcl $eval 15 &
~/Documents/Scripts/Geo_again/Anogeo/Blast_self/4.RBBH_calculate_slave.tcl $eval 16 &
~/Documents/Scripts/Geo_again/Anogeo/Blast_self/4.RBBH_calculate_slave.tcl $eval 17 &
~/Documents/Scripts/Geo_again/Anogeo/Blast_self/4.RBBH_calculate_slave.tcl $eval 18 &
~/Documents/Scripts/Geo_again/Anogeo/Blast_self/4.RBBH_calculate_slave.tcl $eval 19 &
~/Documents/Scripts/Geo_again/Anogeo/Blast_self/4.RBBH_calculate_slave.tcl $eval 20 &

wait
echo "Collating RBBH hits into single Master List ..."

chmod +x ~/Documents/Scripts/Geo_again/Anogeo/Blast_self/4.1.RBBH_collate_slave.tcl
~/Documents/Scripts/Geo_again/Anogeo/Blast_self/4.1.RBBH_collate_slave.tcl $eval &

wait
echo "RBBH: All done"
echo "MCL: Calculating clusters ... "

chmod +x ~/Documents/Scripts/Geo_again/Anogeo/Blast_self/4.2.MCL_slave.sh 
~/Documents/Scripts/Geo_again/Anogeo/Blast_self/4.2.MCL_slave.sh $ival $trunc_eval &

wait
echo "MCL: All done "

