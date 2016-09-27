#!/usr/local/bin/tclsh

#### Download the features tables (if exist), gbb files, and genomic fna files based on the ncbi names of proteome (faa) files in a specific folder ####

source ~/Dropbox/Scripts/General_utils.tcl
set direct /users/aesin/desktop

proc missing_files_announcer {directory} {
	cd $directory
	set out_files [glob -nocomplain *gz]
	if {[llength $out_files] == 0} {puts "No zipped files"; return}
	global num_prok_proteomes
	if {[llength $out_files] == $num_prok_proteomes} {
		puts "\nThe number of files in $directory match the number of prokaryotic proteome files"
	} else {
		set missing_no [expr $num_prok_proteomes - [llength $out_files]]
		puts "\nThere are $missing_no fewer files in $directory than there are prokaryotic proteome files"
	}
}

###########################################################################

set feat_output $direct/Clean_Genomes/Prok_feature_tables
set gbff_output $direct/Clean_Genomes/Prok_gbff_files
set fna_output $direct/Clean_Genomes/Prok_genomic_fna

file mkdir $feat_output
file mkdir $gbff_output
file mkdir $fna_output

cd $direct/Clean_Proteomes/All_fastas_ncbi_id
set proteome_files [glob *faa]

set fail_download_list {}
set num_prok_proteomes 0

foreach proteome_file $proteome_files {

	openfile $direct/Clean_Proteomes/All_fastas_ncbi_id/$proteome_file
	if {[regexp {\(Euk\)} $data] != 1} {

		incr num_prok_proteomes

		regsub {_protein.faa} $proteome_file {} ncbi_id

		set features_file "$ncbi_id\_feature_table.txt.gz"
		regsub {.gz} $features_file {} unz_features_file
		puts "\n--> $features_file"

		set gbff_file "$ncbi_id\_genomic.gbff.gz"
		regsub {.gz} $gbff_file {} unz_gbff_file
		puts "--> $gbff_file"
		
		set genomic_file "$ncbi_id\_genomic.fna.gz"
		regsub {.gz} $genomic_file {} unz_genomic_file
		puts "--> $genomic_file"

		if {[file exists $fna_output/$genomic_file] == 1 || [file exists $fna_output/$unz_genomic_file] == 1} {
			puts "Genomic fna file already downloaded, skipping ..."
			set fna_success 0
		} else {
			cd $fna_output
			puts "Downloading genomic_file ...\n"
			set fna_success [catch {exec /bin/bash -c "curl -O -m 240 -s --retry 10 --retry-delay 1 ftp://ftp.ncbi.nlm.nih.gov/genomes/all/$ncbi_id/$genomic_file"}]
			puts "Fna_success == $fna_success"
		}

		if {$fna_success == 1} {
			puts "Genomic file did not download -> maybe the records is suppresed? Added ncbi_id to non-downloaded list"
			lappend fail_download_list $ncbi_id
			continue
		}

		if {[file exists $gbff_output/$gbff_file] == 1 || [file exists $gbff_output/$unz_gbff_file] == 1} {
			puts "Gbff file already downloaded, skipping ..."
		} else {
			cd $gbff_output
			puts "Downloading gbff_file ..."
			catch {exec /bin/bash -c "curl -O -m 240 -s --retry 10 --retry-delay 1 ftp://ftp.ncbi.nlm.nih.gov/genomes/all/$ncbi_id/$gbff_file"}
		}

		if {[file exists $feat_output/$features_file] == 1 || [file exists $feat_output/$unz_features_file] == 1} {
			puts "Features file already downloaded, skipping ..."
		} else {
			cd $feat_output
			puts "Downloading features_file ..."
			catch {exec /bin/bash -c "curl -O -m 240 -s --retry 10 --retry-delay 1 ftp://ftp.ncbi.nlm.nih.gov/genomes/all/$ncbi_id/$features_file"}
		}

	} else {
		puts "Not prokaryotic - skipping ..."
	}

}
exec clear >@ stdout
puts "######"
missing_files_announcer $fna_output
missing_files_announcer $gbff_output
missing_files_announcer $feat_output
puts "\n######"

puts "\nNumber of completely missing genome files: [llength $fail_download_list]"



