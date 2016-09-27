set direct /scratch/ade110
set org Geo_v_all
set ival 2.0
set eval -100

cd $direct/$org/Rbbh_all_$eval

exec mcl Master_rbh_weight.txt --abc -I $ival -o Master_groups_weight_$ival\.txt -te 32 -scheme 7 --show-log=y --abc-neg-log10 -abc-tf "ceil(200)"

