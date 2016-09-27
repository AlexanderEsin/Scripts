#!/usr/local/bin/tclsh

package require Tk
wm iconify .
###########################################################################

puts "\n\nIs the BLAST all vs all (\"ava\") OR of x proteomes vs y proteomes (\"xvy\")?"
flush stdout
set mode [gets stdin]
regsub -all {"} $mode {} mode

if {$mode != "ava" && $mode != "xvy"} {
	puts "ERROR: You need to select a valid mode\nExiting..."
	exit 2
}

puts "\n...You chose the \"$mode\" mode"

puts "Select or create the directory in which you want the experimental folder to be housed..."
set home_dir [tk_chooseDirectory -initialdir ~/desktop/ -title "Select directory"]
if {$home_dir == ""} {
	puts "\nYou need a directory\nExiting..."
	exit 2
}

if {$mode == "ava"} {
	puts "\nNow choose the directory containing your proteome fasta files..."
	set original_fasta_direct [tk_chooseDirectory -initialdir ~/desktop/ -mustexist TRUE -title "Choose fasta directory"]
	puts "This is your chosen directory ... $original_fasta_direct\n"
	
	file mkdir $home_dir/Input/Fastas
	set input_dir $home_dir/Input
	set fasta_input $home_dir/Input/Fastas

	cd $original_fasta_direct
	set fasta_files [glob -nocomplain *.faa]
	if {[llength $fasta_files] > 0} {
		foreach file $fasta_files {
			file copy -force $file $fasta_input
			puts "$file ... moved to $fasta_input"
		}
	} else {
		puts "You need to have some fasta files!"
	}

	cd $input_dir
	file mkdir DB
	set db_input $input_dir/DB

	cd $fasta_input
	set dblist [glob *faa]
	set i 1

	puts "\n\n"
	foreach db $dblist {
		regsub -all " " $db "_" new_name
		file rename -force $db $new_name
		catch {exec makeblastdb -in $new_name -dbtype prot -out $db_input/$db.db}
		puts "Made database for $new_name ==== $i / [llength $dblist]"
		incr i
	}

	set i [expr $i -1]

	puts "\n===DONE===\nMade databases for: $i sequences\n"
	puts "Your database files are stored in $db_input"
} elseif {$mode == "xvy"} {
	set original_fasta_dirs {}
	puts "\nChoose the FIRST directory containing your proteome (\"x\") fasta files..."
	lappend original_fasta_dirs [tk_chooseDirectory -initialdir ~/desktop/ -mustexist TRUE -title "Choose fasta directory (1)"]
	puts "This is your first chosen directory ... [lindex $original_fasta_dirs 0]\n"

	puts "\nChoose the SECOND directory containing your proteome (\"x\") fasta files..."
	lappend original_fasta_dirs [tk_chooseDirectory -initialdir ~/desktop/ -mustexist TRUE -title "Choose fasta directory (2)"]
	puts "This is your second chosen directory ... [lindex $original_fasta_dirs 1]\n"

	if {[lindex $original_fasta_dirs 0] == [lindex $original_fasta_dirs 1]} {
		puts "BLASTING a directory against itself is equivalent to an ALL v ALL blast - please choose that option\nExiting..."
		exit 2
	}

	set input_dir $home_dir/Input
	file mkdir $input_dir/Fastas_1
	file mkdir $input_dir/Fastas_2
	
	set folder_counter 1
	foreach dir $original_fasta_dirs {
		cd $dir
		set fasta_files [glob -nocomplain *.faa]
		if {[llength $fasta_files] > 0} {
			foreach file $fasta_files {
				file copy -force $file $input_dir/Fastas_$folder_counter
				puts "$file ... moved to $input_dir/Fastas_$folder_counter"
			}
		} else {
			puts "You need to have some fasta files!"
		}
		incr folder_counter
	}

	file mkdir $input_dir/DB_forward
	file mkdir $input_dir/DB_reverse

	cd $input_dir/Fastas_1

	set dblist [glob *faa]
	set i 1
	puts "\n ##### DATABASING FIRST SET #####\n"

	foreach db $dblist {
		regsub -all " " $db "_" new_name
		file rename -force $db $new_name
		catch {exec makeblastdb -in $new_name -dbtype prot -out $input_dir/DB_forward/$db.db}
		puts "Made database for $new_name ==== $i / [llength $dblist]"
		incr i
	}

	set i [expr $i -1]

	cd $input_dir/Fastas_2

	set dblist [glob *faa]
	set x 1
	puts "\n ##### DATABASING SECOND SET #####\n"

	foreach db $dblist {
		regsub -all " " $db "_" new_name
		file rename -force $db $new_name
		catch {exec makeblastdb -in $new_name -dbtype prot -out $input_dir/DB_reverse/$db.db}
		puts "Made database for $new_name ==== $x / [llength $dblist]"
		incr x
	}

	set x [expr $x -1]

	puts "\n===DONE===\nMade databases for Fastas_1: $i sequences\nMade databases for Fastas_2: $x sequences\n"
	puts "Your database files are stored in $input_dir/DB_forward AND $input_dir/DB_reverse"
} else {
	puts "There seems to be an error. Oops. Exiting..."
	exit 2
}

puts "Do you want to run the blast in (p)arallel or on a (s)ingle core?"
flush stdout
set blast_split_switch [gets stdin]


after 1000

exit
