## Rename all the processed genome files to the names that are used in the Species tree so that they match ##

source ~/Dropbox/Scripts/General_utils.tcl

set ncbi_binomial_translation [split [openfile /users/aesin/desktop/Test_ncbi_binomial_transl.tsv] \n]

set direct /users/aesin/desktop/Geo_analysis/Geo_omes/Genome_alignment
cd $direct/Geobac_guide_tree

set tip_list [split [string trim [openfile Tips.txt]] \n]
puts $tip_list

cd $direct/Geobac_genomes_processed
set genomes_files [glob *fna]

puts "Total number of genome files to be renamed: [llength $genomes_files]"

foreach genome_file $genomes_files {
	set ncbi_id [string range $genome_file 0 [string first "genomic" $genome_file]-2]
	puts "\n$ncbi_id"
	set name [lindex [split [lsearch -inline -glob $ncbi_binomial_translation *$ncbi_id] \t] 0]
	regsub -all {\.} $name {} adj_name
	puts $adj_name

	set tip_hit_success [lsearch $tip_list $adj_name]
	if {$tip_hit_success != -1} {
		set tip_hit_match [lsearch -inline $tip_list $adj_name]
		puts "Adjusted name: $adj_name matches a name in the tip list: $tip_hit_match"
		file copy -force $genome_file $direct/Geobac_genomes_renamed/
		file rename -force $direct/Geobac_genomes_renamed/$genome_file $direct/Geobac_genomes_renamed/$tip_hit_match\.fna
	} else {
		foreach tip_name $tip_list {
			set reverse_match [lsearch -all -inline -glob $adj_name $tip_name*]
			if {[llength $reverse_match] == 1} {
				puts "The corresponding tip name in the species tree does not match exactly, BUT the corresponding name was found as: $tip_name\nDoes this correspond to the adjusted name: $adj_name\?"
				file copy -force $genome_file $direct/Geobac_genomes_renamed/
				file rename -force $direct/Geobac_genomes_renamed/$genome_file $direct/Geobac_genomes_renamed/$tip_name\.fna
			}
		}
	}
}
