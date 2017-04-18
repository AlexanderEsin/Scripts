#!/usr/local/bin/tclsh
### Split equally, compress, and upload the query files for the DB vs QUERY blast direction ###

###########################################################################
set direct /users/aesin/desktop/Consensus_trees
###
set org Inside_group
set input_dir $direct/$org/Input
set query_dir Fastas
#n=number of cores wanted - 1 (e.g. for 64 cores n = 63)
set n 31
### Zip(1) or no zip(0)?
set zip 1
### Scp(1) or no scp(0)?
set scp 1
###########################################################################

puts "\nSplitting into $n directories ..."
cd $input_dir

file mkdir Split

cd $input_dir/$query_dir
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
		file mkdir $input_dir/Split/$dir
		puts "Directory made: $dir out of [expr $n + 1]"
		set counter $i
		while {$counter < [llength $org_size_list]} {
			file copy [lindex $org_size_list $counter] $input_dir/Split/$dir/
			set counter [expr ($counter + ($n + 1))]
			incr y
		}
		incr i
	}
}

puts "\nSplit done"

###########################################################################

if {$zip == 1} {
	cd $direct/$org
	puts "\nNow zipping ..."
	exec pigz -q -r --fast Input
}
if {$scp == 1} {
	puts "\nNow uploading to ax3 ..."
	cd $direct/$org
	exec rsync -r Input ade110@ax3.hpc.ic.ac.uk:/csc/rawdata/Warnecke/ADE/Consensus_trees/$org >&@stdout
	# This allows us to execute the gunzip command on the server. Note the use of the quotation marks around the unix command - since single quotes have no special meaning in tcl
	exec ssh ade110@ax3.hpc.ic.ac.uk "gunzip -r /csc/rawdata/Warnecke/ADE/Consensus_trees/$org/Input"
}

puts "\nUpload to server and unzip done"
###########################################################################
