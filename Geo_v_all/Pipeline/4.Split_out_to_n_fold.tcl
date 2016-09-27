### This script splits the accumulated blast hit files (in out_all) into a number of folders spcecified by (n-1). The blast hits are paired for downstream Rbbh analysis. ###

set direct /scratch/ade110
set org Geo_v_all

#n=number of split folders wanted - 1
set n 31
###

###########################################################################

cd $direct/$org

set i 0

file mkdir Out_split
cd Out_split
while {$i <= $n} {
	set dir [expr ($i+1)]
	file mkdir $dir
	incr i
}

cd $direct/$org/Out_all
set orig_list [glob *.tsv]
set blasthit $orig_list

set i 0
set y 0
set file_counter 0

while {[llength $blasthit] > $y} {
	if {$i == [expr ($n+1)]} {
		set i 0
	} else {
		set dir [expr ($i+1)]

		set matchA [lindex $blasthit $file_counter]
		if {[file exists $matchA] == 1} {
			set rearrange [split $matchA &]
			set q [lindex $rearrange 0]
			set s [lindex $rearrange 1]
			regsub .tsv $s {} s
			set matchB "$s&$q\.tsv"

			file rename $matchA $direct/$org/Out_split/$dir
			incr y
			file rename $matchB $direct/$org/Out_split/$dir
			incr y
			set file_counter [expr $file_counter + 1]
			puts "$y / [llength $orig_list]"
			incr i
		} else {
			set file_counter [expr $file_counter + 1]
		}
	}
}

puts "\n\n$y"

puts "====== DONE ======"

# cd /scratch/ade110/Geo_v_all/Out_split
# set i 1
# set dirs [glob -type d *]
# foreach dir $dirs {
# 	cd /scratch/ade110/Geo_v_all/Out_split/$dir
# 	set files [glob -nocomplain *.tsv]
# 	foreach file $files {
# 		file rename $file /scratch/ade110/Geo_v_all/Out_all
# 		puts "$i"
# 		incr i
# 	}
# }

# # set total_length 0
# # set dirs [glob -type d *]
# # foreach dir $dirs {
# # 	cd /scratch/ade110/Consensus_trees/Bacillaceae/Out_split/$dir
# # 	set files [glob -nocomplain *.tsv]
# # 	set total_length [expr $total_length + [llength $files]]
# # }
# # puts "NO: $total_length"
