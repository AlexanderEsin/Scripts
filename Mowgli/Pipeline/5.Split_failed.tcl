#!/usr/bin/tclsh
## Split the remaining (Failed) files into a n folders to serve as input for the 64-bit Mac Pro Mowgli ##

set transfer_cost 3
set input_type folders; # "files" or "folders"
set file_ext ".txt"

###########################################################################

### Assign variables ###
set direct /users/aesin/desktop/Mowgli
set input_dir_name MPRO_input/Input_BS_$transfer_cost
set output_dir_name MPRO_input/MPRO_input_BS_$transfer_cost

# n = number of folders to split into #
set n 2

###########################################################################

### Prepare the unsplit files and directories ###
cd $direct
file mkdir $output_dir_name
cd $direct/$input_dir_name

if {$input_type eq "files"} {

	set unsplit_files [glob *$file_ext]
	
	### Set the counters ###
	set i 0
	set y 1

	### Do the splitting ###
	while {$y <= [llength $unsplit_files]} {
		while {$i <= [expr $n - 1]} {
			set dir [expr ($i + 1)]
			file mkdir $direct/$output_dir_name/$dir
			puts "Directory made: $dir out of [expr $n]"
			set counter $i
			while {$counter < [llength $unsplit_files]} {
				file copy [lindex $unsplit_files $counter] $direct/$output_dir_name/$dir
				set counter [expr ($counter + $n)]
				incr y
			}
			incr i
		}
	}

	puts "All directories made"
	exit

} elseif {$input_type eq "folders"} {

	set unsplit_folders [glob -type d *]

	### Set the counters ###
	set i 0
	set y 1

	### Do the splitting ###
	while {$y <= [llength $unsplit_folders]} {
		while {$i <= [expr $n - 1]} {
			set dir [expr ($i + 1)]
			file mkdir $direct/$output_dir_name/$dir
			puts "Directory made: $dir out of [expr $n]"
			set counter $i
			while {$counter < [llength $unsplit_folders]} {
				cd $direct/$input_dir_name/[lindex $unsplit_folders $counter]
				set file [lindex [glob *$file_ext] 0]
				file copy $file $direct/$output_dir_name/$dir
				set counter [expr ($counter + $n)]
				incr y
			}
			incr i
		}
	}

	puts "All directories made"
	exit

}

