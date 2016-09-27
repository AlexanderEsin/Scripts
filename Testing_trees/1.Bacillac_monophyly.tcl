set direct /users/aesin/desktop
set ival 2.0

proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}

package require sqlite3
sqlite3 tax_ranks $direct/Clean_proteomes/Taxonomic_ranks_db

### Set required taxonomic groups ###
set umbrella_taxon_group Firmicutes
set subset_taxon_groups {Bacilli Bacillaceae Bacillales}	
###########################################################################

### SWITCHES ###
set folder_make_switch 1
set rscript_only_switch 0

### FOLDER MAKING ####
if {$folder_make_switch == 1} {
	file mkdir /users/aesin/desktop/Tree_sorting/Tree_build
	file mkdir /users/aesin/desktop/Tree_sorting/Final_trees/Final_BS
	file mkdir /users/aesin/desktop/Tree_sorting/Final_trees/Analysis
	file mkdir /users/aesin/desktop/Tree_sorting/Unfinished_trees/
	

	cd $final_tree_dir
	file mkdir Bacillac_mono/All_bacillac_anoxygeo
	file mkdir Bacillac_mono/Not_all_bacillac_anoxygeo/All_bacillac_bacillales
	file mkdir Bacillac_mono/Not_all_bacillac_anoxygeo/Not_all_bacillac_bacillales

	file mkdir Bacillac_not_mono/Geo_Anoxy_mono/Bacillales_mono
	file mkdir Bacillac_not_mono/Geo_Anoxy_mono/Bacillales_not_mono
	file mkdir Bacillac_not_mono/Geo_Anoxy_not_mono

	file mkdir Bacillac_only/Geo_Anoxy_only
	file mkdir Bacillac_only/Geo_Anoxy_mono
	file mkdir Bacillac_only/Geo_Anoxy_not_mono

	file mkdir KEYS
}

set geo_v_all_dir $direct/Geo_v_all/$ival
set tree_build_dir $direct/Tree_sorting/Tree_build
set bs_final_dir $direct/Tree_sorting/Final_trees/Final_BS
set final_tree_dir $direct/Tree_sorting/Final_trees/Analysis
set unfinished_dir $direct/Tree_sorting/Unfinished_trees

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

set folders [glob -type d *]
if {$rscript_only_switch == 0} {
	foreach folder $folders {
		cd $tree_build_dir/$folder
		puts "Making BS-stripped ML tree and pulling out the requested tip files ... $i / [llength $folders]"

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
			continue
		} else {
			file copy -force $folder\_tree.txt $bs_final_dir
		}

		set besttrees [glob -nocomplain *bestTree.*]
		if {[llength $besttrees] > 0} {
			foreach tree $besttrees {
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
		}

		if {[llength $umbrella_taxon_group] != 1} {
			puts "THERE CAN BE ONLY ONE! ... umbrella taxon group.. -> fix yo life"
			exit 2
		}
		
		if {[file exists $folder\_$umbrella_taxon_group\_tips.txt] == 0} {
			set umbrella_tips {}
			set umbrella_entries {}

			set key_entries [split [string trim $key] \n]
			foreach entry $key_entries {
				set binomial [lindex [split $entry \t] 1]

				tax_ranks eval {SELECT ranks FROM t1 WHERE binomial = $binomial} {
					set split_ranks [split $ranks \t]
				}

				if {[lsearch $split_ranks $umbrella_taxon_group] > -1} {
					lappend umbrella_tips "\"[lindex $entry 1]_\{[lindex $entry 2]\}_\{[lindex $entry 3]\}\""
					lappend umbrella_entries $entry
				}

			}
			set umbrella_tips [join $umbrella_tips "\n"]
			set out [open $folder\_$umbrella_taxon_group\_tips.txt w]
			puts $out [string trim $umbrella_tips]
			close $out
		} else {
			set umbrella_entries [split [string trim $key] \n]
		}
		
		foreach group $subset_taxon_groups {
			if {[file exists $folder\_$group\_tips.txt] == 0} {
				set subset_tips {}
				foreach entry $umbrella_entries {
					set binomial [lindex [split $entry \t] 1]

					tax_ranks eval {SELECT ranks FROM t1 WHERE binomial = $binomial} {
						set split_ranks [split $ranks \t]
					}

					if {[lsearch $split_ranks $group] > -1} {
						lappend subset_tips "\"[lindex $entry 1]_\{[lindex $entry 2]\}_\{[lindex $entry 3]\}\""
					}
				}
				set subset_tips [join $subset_tips "\n"]
				set out [open $folder\_$group\_tips.txt w]
				puts $out [string trim $subset_tips]
				close $out
			}
		}
		incr i
	}
}

tax_ranks close

# Open up the R script template and substitute the correct working directories #
openfile /users/aesin/desktop/Scripts/Testing_trees/0.MProot_tree_bacillac.R
regsub {TREE_BUILD} $data $tree_build_dir temp_data
regsub {FINAL_TREE} $temp_data $final_tree_dir temp_data
#puts $temp_data

set out [open /users/aesin/desktop/scripts/Testing_trees/0.MProot_tree_bacillac_temp.R w]
puts $out $temp_data
close $out

# Run the R script - it should pop up with a little Tk progress bar (fucking around) #
file attributes /users/aesin/desktop/scripts/Testing_trees/0.MProot_tree_bacillac_temp.R -permissions 0777
## THIS REDIRECT TO STDOUT IS THE FUCKING HOLY GRAIL. YOU SPENT 2 HOURS TRYING TO FIGURE OUT HOW TO DO THIS ##
exec /users/aesin/desktop/scripts/Testing_trees/0.MProot_tree_bacillac_temp.R >@stdout
#file delete /users/aesin/desktop/scripts/Testing_trees/0.MProot_tree_bacillac_temp.R

puts "\n\n==== DONE ===="


### CLEAN TROUBLE SHOOTING ###

# set tree_build_dir /users/aesin/desktop/Tree_sorting/Tree_build_2.5

# proc openfile {fl} {
# 	global data
# 	set in [open $fl r]
# 	set data [read $in]
# 	close $in
# 	return
# }

###########################################################################

# cd $tree_build_dir

# set folders [glob -type d *]

# foreach folder $folders {
# 	cd $tree_build_dir/$folder
# 	set tip_files [glob -nocomplain *tips.txt]
# 	set ana_file [lindex [glob -nocomplain *ana.txt*] 0]
# 	set mp_file [lindex [glob -nocomplain *mp_tree.txt*] 0]
# 	file delete $ana_file
# 	file delete $mp_file
# 	foreach file $tip_files {
# 		file delete $file
# 	}
# }
