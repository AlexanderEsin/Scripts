
#Change the directory below
set direct /scratch/ade110
set org Geo_v_all
###

proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}

#Refill the original out directory with file from temp
cd $direct/$org/Out_split/[lindex $argv 0]/temp
set files [glob *tsv]
foreach file $files {
	file rename $file $direct/$org/Out_split/[lindex $argv 0]/$file
}
cd ..
file delete -force temp

