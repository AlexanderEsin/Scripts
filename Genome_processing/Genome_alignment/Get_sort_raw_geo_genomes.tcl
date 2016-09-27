#!/usr/local/bin/tclsh

# ### Get all the Geobacillus genomes from the genome set downloaded from NCBI (these may be completed genomes or contig assemblies) and transfer them into a directory in Geo_analysis ###
source ~/Dropbox/Scripts/General_utils.tcl

set ncbi_binomial_translation [split [openfile /users/aesin/desktop/Test_ncbi_binomial_transl.tsv] \n]

# set geo_counter 0
# foreach entry $ncbi_binomial_translation {
# 	set name_tag [lindex [split $entry \t] 0]
# 	if {[string match "Geobacillus*" $entry] == 1} {
# 		incr geo_counter
# 		set ncbi_id [lindex [split $entry \t] 1]
# 		set genomic_fna_file_name "$ncbi_id\_genomic.fna"
# 		set genomic_gbff_file_name "$ncbi_id\_genomic.gbff"
# 		## Copy fna files ##
# 		#file copy -force /users/aesin/desktop/Clean_genomes/Prok_genomic_fna/$genomic_fna_file_name /users/aesin/desktop/Geo_analysis/Geo_omes/Geobac_genomes_raw
# 		## Cope gbff files ##
# 		file copy -force /users/aesin/desktop/Clean_genomes/Prok_gbff_files/$genomic_gbff_file_name /users/aesin/desktop/Geo_analysis/Geo_omes/Geobac_genomes_genbank
# 	}
# }

# puts "\n\n$geo_counter"

### Split the genomes into those that are fully assembled versus those that are contig assemblies. Those that are fully assembled (n genetic elements <= 5) are search for elements with the word plasmid, and the plasmid sequence is removed. ###

set direct /users/aesin/desktop/Geo_analysis/Geo_omes/Genome_processing
file mkdir $direct/Geobac_genomes_assembled_no_plasmids
file mkdir $direct/Geobac_genomes_contigs
cd $direct/Geobac_genomes_raw
set genomes [glob *fna]

foreach genome $genomes {
	openfile $genome
	set genomic_elements [regexp -all ">" $data]
	puts "$genome\t\t$genomic_elements"

	regsub "_genomic.fna" $genome {} ncbi_id
	set name [lindex [split [lsearch -inline -glob $ncbi_binomial_translation *$ncbi_id] \t] 0]


	if {$genomic_elements <= 5} {
		set no_plasmid_genome {}
		set elements [split_genes $data]
		foreach element $elements {
			set element_header [string range $element 0 [string first \n $element]]
			#puts $element_header
			if {[regexp "plasmid" $element_header] != 1} {
				lappend no_plasmid_genome $element
			}
		}

		set no_plasmid_genome [join $no_plasmid_genome]

		set no_plasmid_genome_elements [regexp -all ">" $no_plasmid_genome]
		#puts "$genome\t$no_plasmid_genome_elements"

		if {$no_plasmid_genome_elements == 1} {
			#puts "Fully assembled: $genome\t$name"
			# set out [open $direct/Geobac_genomes_assembled_no_plasmids/$genome w]
			# puts $out $no_plasmid_genome
			# close $out
		} else {
			#puts "Plasmid not identified: $name\t$genome\t$no_plasmid_genome_elements"
		}
	} else {
		#file copy -force $genome $direct/Geobac_genomes_contigs
		#puts "Contig genome: $genome\t$name\t$genomic_elements"
	}

}
