#!/usr/local/bin/tclsh
set taxon [lindex $argv 0]
set direct /users/aesin/desktop
set org Consensus_trees/Rooting/$taxon
set ival 5
set eval -10

cd $direct/$org/Rbbh_all_$eval

exec mcl Master_rbh_weight.txt --abc -I $ival -o Master_groups_weight_$ival\.txt -te 8 -scheme 7 --show-log=y --abc-neg-log10 -abc-tf "ceil(200)" >&@stdout
