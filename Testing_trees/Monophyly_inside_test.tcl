#!/usr/local/bin/tclsh
set direct /users/aesin/desktop
set inside_group_folder $direct/Clean_Proteomes/Inside_group
set ival 2.0
set desktop_dir Tree_sorting

###########################################################################

source ~/desktop/Scripts/General_utils.tcl

###########################################################################

### SWITCHES ###
set folder_make_switch 1
set make_tips_mp_tree_switch 1

### FOLDER MAKING ####
set geo_v_all_dir $direct/Geo_v_all/$ival
set tree_build_dir $direct/$desktop_dir/Tree_build
set bs_final_dir $direct/$desktop_dir/Final_trees/Final_BS
set final_tree_dir $direct/$desktop_dir/Final_trees/Analysis
set unfinished_dir $direct/$desktop_dir/Unfinished_trees

if {$folder_make_switch == 1} {
	file mkdir $tree_build_dir
	file mkdir $bs_final_dir
	file mkdir $final_tree_dir 
	file mkdir $unfinished_dir
	

	cd $final_tree_dir
	file mkdir Anoxygeo_mono/Anoxy_geo_only
	file mkdir Anoxygeo_mono/No_inside_group
	file mkdir Anoxygeo_mono/All_inside_group
	file mkdir Anoxygeo_mono/Inside_group
	file mkdir Anoxygeo_not_mono/Geo_mono
	file mkdir Anoxygeo_not_mono/Geo_not_mono

	file mkdir KEYS
}

### Copy the Tree_build directory from Geo_v_all to the $tree_build_dir ###
cd $geo_v_all_dir/Tree_build
set sizes [glob -type d *]
set i 0
foreach size $sizes {
	cd $geo_v_all_dir/Tree_build/$size
	set folders [glob -nocomplain -type d *]
	foreach folder $folders {
		if {[file exists $tree_build_dir/$folder] == 0} {
			file copy -force $folder $tree_build_dir/
			file copy -force $geo_v_all_dir/Group_fastas_MSA/KEY_$folder\.txt $tree_build_dir/$folder/
			incr i
			puts "Copying folders .... $i"
		}
	}
}

cd $tree_build_dir
set i 1

set pata {k+?[0-9]+?:+?}
set patb {\t+?.+?\n+?}
set patc {\t+.+?}

set folders [lsort -dictionary [glob -type d *]]
if {$make_tips_mp_tree_switch == 1} {

	cd $inside_group_folder
	set inside_group_proteomes [glob *.faa]
	openfile $direct/Clean_Proteomes/All_fastas_geo_mark_key_file.tsv
	set geo_mark_keys [split [string trim $data] \n]
	foreach inside_group $inside_group_proteomes {
		set inside_group_hit [lsearch -inline -glob $geo_mark_keys $inside_group*]
		set inside_group_id [lindex [split $inside_group_hit \t] end]
		lappend inside_group_id_l $inside_group_id
	}

	set failed_tree_build_l {}

	foreach folder $folders {
		cd $tree_build_dir/$folder
		puts "Making BS-stripped ML tree and pulling out the necessary tip files ... $i / [llength $folders]"

		###
		# In the tcl portion, we prep a tree file with the correct tip labels, but not including the BS values, for downstream monophyly and branching calculations. It is essentially identical to the process at the back end of the RAxML scripts #
		set key [lindex [glob KEY*] 0]
		file copy -force $key $final_tree_dir/KEYS
		openfile $key
		set key $data

		if {[file exists $folder\_tree.txt] == 0} {
			puts "\tERROR: the tree file for folder $folder does not exist here. Moving the folder to $unfinished_dir ..."
			cd ..
			if {[file exists $unfinished_dir/$folder] == 1} {
				file delete -force $tree_build_dir/$folder
			} else {
				file rename -force $tree_build_dir/$folder $unfinished_dir/
			}
			incr i
			lappend failed_tree_build_l $folder
			continue
		} else {
			file copy -force $folder\_tree.txt $bs_final_dir
		}

		set besttrees [glob -nocomplain *bestTree.*]
		set tree [lindex $besttrees 0]

		if {[llength $tree] > 0} {
			openfile $tree
			set branchids [regexp -all -inline $pata $data]
			foreach k_val $branchids {
				regsub {:} $k_val {} k_val_trunc
				regexp -line $k_val_trunc$patb $key hit
				set binomial [lindex [split [string trim $hit] \t] 1]
				set class [lindex [split [string trim $hit] \t] 2]
				set id [lindex [split [string trim $hit] \t] 3]
				# If you use k_val_trunc below, you may substitute k17 when you want to find k1 - so searching for k1: is more accurate #
				regsub $k_val $data "$binomial\_\{$class\}_\{$id\}:" data
			}
			set out [open $folder\_mp_tree.txt w]
			puts $out $data
			close $out
		}


		if {[file exists $folder\_inside_group_tips.txt] == 0} {
			set inside_group_tips_l {}

			set key_entries [split [string trim $key] \n]
			foreach entry $key_entries {
				set id [lindex [split $entry \t] end]
				set species_id [string range $id [string last "." $id]+1 end]

				if {[lsearch $inside_group_id_l $species_id] > -1} {
					lappend inside_group_tips_l "\"[lindex $entry 1]_\{[lindex $entry 2]\}_\{[lindex $entry 3]\}\""
				}
			}

			set inside_group_tips_str [join $inside_group_tips_l "\n"]
			set out [open $folder\_inside_group_tips.txt w]
			puts $out [string trim $inside_group_tips_str]
			close $out
		}
		
		incr i
	}
}

# Open up the R script template and substitute the correct working directories #
openfile /users/aesin/desktop/Scripts/Testing_trees/Inside_test_R.R
regsub {TREE_BUILD} $data $tree_build_dir temp_data
regsub {FINAL_TREE} $temp_data $final_tree_dir temp_data
#puts $temp_data

set out [open /users/aesin/desktop/scripts/Testing_trees/Inside_test_R_temp.R w]
puts $out $temp_data
close $out

# # Run the R script #
file attributes /users/aesin/desktop/scripts/Testing_trees/Inside_test_R_temp.R -permissions 0777
# ## THIS REDIRECT TO STDOUT IS THE FUCKING HOLY GRAIL. YOU SPENT 2 HOURS TRYING TO FIGURE OUT HOW TO DO THIS ##
exec /users/aesin/desktop/scripts/Testing_trees/Inside_test_R_temp.R >@stdout
file delete /users/aesin/desktop/scripts/Testing_trees/Inside_test_R_temp.R

puts "\n\n==== DONE ===="
if {$rscript_only_switch == 0} {
	set failed_tree_build_str [join $failed_tree_build_l \n]
	puts "There were: [llength $failed_tree_build_l] incomplete tree reconstructions. These were:\n$failed_tree_build_str"
}
