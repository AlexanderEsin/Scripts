### Split equally, compress, and upload the query files for the DB vs QUERY blast direction ###

###########################################################################
set direct /users/aesin/desktop
###
set org Geo_v_all/Input
set query_dir All_fastas
#n=number of cores wanted - 1 (e.g. for 64 cores n = 63)
set n 127
###########################################################################

cd $direct/$org

file mkdir Split

cd $direct/$org/$query_dir
set unsortlist [exec ls -s -S]
set header [string first \n $unsortlist]
set unsortlist [string trim [string range $unsortlist $header end]]

set pata {.+?[0-9]\s+?}
regsub -all -line $pata $unsortlist {} sortlist
set org_size_list [string trim $sortlist]

set i 0
set y 1

while {$y <= [llength $org_size_list]} {
	while {$i <= $n} {
		set dir [expr ($i + 1)]
		file mkdir $direct/$org/Split/$dir
		puts "Directory made: $dir out of [expr $n + 1]"
		set counter $i
		while {$counter < [llength $org_size_list]} {
			file copy [lindex $org_size_list $counter] $direct/$org/Split/$dir/
			set counter [expr ($counter + ($n + 1))]
			incr y
		}
		incr i
	}
}

###########################################################################

cd $direct/Geo_v_all
puts "\nNow zipping"
pigz -q -r --fast Input
puts "DONE"

#scp -r Split ade110@ax3.hpc.ic.ac.uk/scratch/ade110/
