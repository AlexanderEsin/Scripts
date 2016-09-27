set direct /scratch/ade110/Geo_v_all/TT
set tree_no 160
set core_no 4

###########################################################################

cd $direct/$tree_no

set fas [lindex [glob *fas] 0]

regsub {.fas} $align {} out_number

exec raxml -f a -s $align -n $out_number\.txt -m PROTCATAUTO -p 1234 -x 10001 -N 100 -T $core_no



