set direct /csc/rawdata/Warnecke/ADE
set org Consensus_trees/Inside_group

set ival 5
set eval -10

cd $direct/$org/Rbbh_all_$eval

exec mcl Master_rbh_weight.txt --abc -I $ival -o Master_groups_weight_$ival\.txt -te 4 -scheme 7 --show-log=y --abc-neg-log10 -abc-tf "ceil(200)" >& Log_MCL_$ival\.txt

