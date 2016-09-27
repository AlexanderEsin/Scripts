set direct /users/aesin/desktop/Consensus_trees
set eval -10
set ival 2
set group Inside_group/$eval

## Raxml options
set bootstrap_no 100
set thread_no 7
set model_name PROTCATAUTO

proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}

###########################################################################

cd $direct/$group
file mkdir Individual_trees_$ival
file mkdir Final_trees_$ival/Individual_trees

cd $direct/$group/Group_fastas_MSA_2
set aligns [glob *.fas]

## Create a new folder for each of the alignment in the Individual trees directory, copy the Species key into each one ##
foreach alignment $aligns {
	set group_number [string range $align 0 end-4]
	file mkdir $direct/$group/Individual_trees_$ival/$group_number
	file copy $alignment $direct/$group/Individual_trees_$ival/$group_number
	file copy Species_KEY.txt $direct/$group/Individual_trees_$ival/$group_number
}


cd $direct/$group/Individual_trees_$ival
set alignment_dirs [glob -type d *]
puts "Will now process [llength $alignment_dirs] individual alignments..."

foreach dir $alignment_dirs {
	cd $direct/$group/Individual_trees_$ival/$dir
	set align_file [lindex [glob *fas] 0]
	set filename [string range $align_file 0 end-4]
	
	## Not sure what this block does ##
	# openfile $align_file
	# set new_list [string trim $data]
	# regsub -all "\n\n" $new_list "\n" new_list
	# set out [open $align w]
	# puts $out $new_list
	# close $out
	# regsub {.fas} $align {} filename
	# file mkdir $direct/Trees/$filename
	# file copy Species_KEY.tsv $direct/Trees/$filename/

	##############################
	puts "---> Now working on tree $filename"
	catch {exec raxml -f a -s $align_file -n $filename.txt -m $model_name -p 1234 -x 10001 -N $bootstrap_no -T $thread_no}

	## Delete the reduced alignment file if it exists ##
	set reduced_file [lindex [glob -nocomplain *.reduced] 0]
	if {[llength $reduced_file] == 1} {
		file delete $reduced_file
	}

	## Get the key information ##
	set key_file [glob *KEY*]
	openfile $key; set key $data

	set pata {k+?[0-9]+?:+?}
	set patb {\t+.+}

	set besttrees [glob -nocomplain *bipartitions.*]
	if {[llength $besttrees] > 0} {
		foreach tree $besttrees {
			openfile $tree
			set branchids [regexp -all -inline $pata $data]
			foreach k_val $branchids {
				regsub {:} $k_val {} k_val_trunc
				regexp -line $k_val_trunc$patb $key hit
				set binomial [lindex [split [string trim $hit] \t] 1]
				# set id [lindex [split [string trim $hit] \t] 2]
				# If you use k_val_trunc below, you may substitute k17 when you want to find k1 - so searching for k1: is more accurate #
				regsub $k_val $data "$binomial\:" data
			}
			set out [open $filename\_tree.txt w]
			puts $out $data
			close $out

			file copy -force $filename\_tree.txt $direct/$group/Final_trees_$ival/Individual_trees
		}
	} else {
		puts "ERROR: Bipartitions file not found"
		exit 2
	}
}

















