set direct /users/aesin/desktop
set org Geobac/Outgroups/1_out_subtilis

proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}

###########################################################################

cd $direct/$org/Out_all

set files [glob *.tsv]
foreach file $files {
	if {[regexp {Geobacillus_sp._A8} $file] == 1} {
		file rename $file /scratch/ade110/Geo_v_all/Out_nw/
	} elseif {[regexp {Geobacillus_sp._CAMR5420} $file] == 1} {
		file rename $file /scratch/ade110/Geo_v_all/Out_nw/
	} elseif {[regexp {Geobacillus_sp._FW23} $file] == 1} {
		file rename $file /scratch/ade110/Geo_v_all/Out_nw/
	}
}