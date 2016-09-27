#!/usr/local/bin/tclsh

### Split equally, compress, and upload the query files for the DB vs QUERY blast direction ###

###########################################################################
set direct /users/aesin/desktop/Consensus_trees/Rooting
###
set org [lindex $argv 0]
set input_dir $org/Input
puts $input_dir
set query_dir DB_forward
#n=number of cores wanted - 1 (e.g. for 64 cores n = 63)
set n 7
###########################################################################

cd $direct/$input_dir

file mkdir DB_split

cd $direct/$input_dir/$query_dir
set sortlist [glob *db.pin]

set org_size_list [string trim $sortlist]

set i 0
set y 1

while {$y <= [llength $org_size_list]} {
	while {$i <= $n} {
		set dir [expr ($i + 1)]
		file mkdir $direct/$input_dir/DB_split/$dir
		puts "Directory made: $dir out of [expr $n + 1]"
		set counter $i
		while {$counter < [llength $org_size_list]} {
			set db_file [lindex $org_size_list $counter]
			regsub {\.pin} $db_file {} db_file
			set db_files [glob $db_file\*]
			foreach file $db_files {
				file copy $file $direct/$input_dir/DB_split/$dir/
			}
			set counter [expr ($counter + ($n + 1))]
			incr y
		}
		incr i
	}
}

###########################################################################

exec /users/aesin/desktop/scripts/Consensus_trees/Rooting/3.HPC.Blastp_master.sh $org

###########################################################################

# cd $direct/$org
# cd ..
# puts "\nNow zipping"
# exec pigz -q -r --fast Input
# puts "DONE"

#scp -r Split ade110@ax3.hpc.ic.ac.uk/scratch/ade110/
