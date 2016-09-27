set ival 2.0
set direct /scratch/ade110/Geo_v_all/$ival

proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}

###########################################################################\

openfile $direct/Distinct_align_patterns.tsv
set align_pat_info [split [string trim $data] \n]

cd $direct/Tree_build

set sizes [glob -type d *]
foreach size $sizes {
	set i 1
	cd $direct/Tree_build/$size
	set tree_dirs [glob -nocomplain -type d *]
	foreach dir $tree_dirs {
		cd $direct/Tree_build/$size/$dir
		# Check whether the final tree file exists (with the tip names having already been exchanged) #
		if {[file exists $dir\_tree.txt] == 0} {

			#puts $direct/Tree_build/$size/$dir

			# Remove all the excess RAxML files #
			set tree_build_files [glob -nocomplain RAxML*]
			foreach file $tree_build_files {
				file delete $file
			}

			# Use the alignment_pattern table to figure out to which core_split this file belongs #
			set align_pat_entry [lsearch -all -inline -glob $align_pat_info $dir\t*]
			if {[llength $align_pat_entry] != 1} {puts "Woof error"; exit 2}
			set pat_number [lindex [split [lindex $align_pat_entry 0] \t] 1]
			if {$pat_number > 4000} {
				set core_split individual
			} elseif {$pat_number > 1250} {
				set core_split 8_core
			} else {
				set core_split 4_core
			}

			# Within the correct core_split and size, pick the smallest split folder #
			cd $direct/MSA_split/$core_split/$size
			set input_folders [glob -type d *]
			set smallest_folder ""
			set smallest_folder_len 1000
			foreach folder $input_folders {
				cd $folder
				set len_folder [llength [glob -nocomplain -type d *]]
				if {$len_folder < $smallest_folder_len} {
					set smallest_folder $folder
					set smallest_folder_len $len_folder
				}
				cd ..
			}

			# Send the unfinished tree build folder to the smallest folder #
			file rename $direct/Tree_build/$size/$dir $direct/MSA_split/$core_split/$size/$smallest_folder/
			puts "$size folder: $dir was moved back to $direct/MSA_split/$core_split/$size/$smallest_folder"

			cd $direct/Tree_build/$size
		} else {
			cd ..
		}
		puts "$size ==== $i / [llength $tree_dirs]"
		incr i
	}
}