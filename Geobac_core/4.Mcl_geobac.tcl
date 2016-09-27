set direct /users/aesin/desktop/Consensus_trees/Geobac/All_geobac
set ival 5
set eval -40

cd $direct
set Rbbh_dirs [glob Rbbh_all*]

foreach dir $Rbbh_dirs {
	cd $direct/$dir
	if {[file exists RBBH_Hitfile_$eval\.txt] == 1} {
		mcl Master_rbh_weight.txt --abc -I $ival -o Master_groups_weight_$ival\.txt -te 7 -scheme 7 --show-log=y --abc-neg-log10 -abc-tf "ceil(200)"
	}
}



######

#-te corresponds to no. threads used
#-scheme (the higher the number the more refined the analysis)
#-abc-neg-log10 allows to use log10 values as edge weights
#-abc-tf (ceil200) caps the highest weight at 200 when 1E-199 = 199