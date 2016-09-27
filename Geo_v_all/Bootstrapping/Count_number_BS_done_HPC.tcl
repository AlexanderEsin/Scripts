source ~/Scripts/General_utils.tcl

set direct /scratch2/ade110/Geo_v_all/2.0/Trees/Final_NBS

cd $direct
set folders [glob -type d *]
puts "There are [llength $folders] trees for BS reconstruction"

set all_completed 0
foreach folder $folders {

	cd $direct/$folder/Bootstraps

	set bs_files [glob -nocomplain BS*]
	if {[llength $bs_files] < 10} {
		puts "There are fewer than 10 BS files in folder: $folder. There are only [llength $bs_files] files."
		puts "\nThe following files are completed:\n[join $bs_files \n]\n"
		continue
	}

	set master_BS_file {}
	foreach bs_file $bs_files {
		set 10_bootstraps [string trim [openfile $bs_file]]
		set file_bs_counter [regexp -all ";" $10_bootstraps]
		if {$file_bs_counter < 10} {
			puts "Folder $folder\: In BS file $bs_file there were only $file_bs_counter bootstraps."
			continue
		}
	}

	incr all_completed
}

puts "$all_completed out of [llength $folders] have completed 100 BS replicates."