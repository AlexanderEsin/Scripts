#!/usr/local/bin/tclsh
set direct /users/aesin/desktop
set org Consensus_trees/Geobac/Outgroups/With_excl_geos

set ival 2
set eval -10

cd $direct/$org/Rbbh_all_$eval

exec mcl Master_rbh_weight.txt --abc -I $ival -o Master_groups_weight_$ival\.txt -te 4 -scheme 7 --show-log=y --abc-neg-log10 -abc-tf "ceil(200)" >& Log_MCL_$ival\.txt

