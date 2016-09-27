### Assign variables ###
set direct /users/aesin/desktop/Deer/Assembly
set input_dir_name Alpha_input
set output_dir_name Alpha_input_split
set file_ext ".fasta"
# n = number of folders to split into #
set n 4

### Prepare the unsplite files and directories ###
cd $direct
file mkdir $output_dir_name
cd $direct/$input_dir_name
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