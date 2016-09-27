#!/usr/local/bin/tclsh
######
# Remember to change the file type to be downloaded in the shell script #
#####
set org_type "Bacteria_Curated"; #FOLDER NAME / Name of organism group e.g. Archaea. This must match a REF file.
set filetype "genomic.gbff.gz"; #Filetype suffix needed
set direct /users/aesin/desktop
######


proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}

###########################################################################
#Copy the assembly_summary_refseq file into the new directory. Also copy the relevant reference list from the Refs folder
#The ref file is just a recent manually copied list of all the relevant species names from the ftp site

file mkdir $direct/$org_type
file copy -force $direct/Scripts/Gen_Pull/assembly_summary_refseq.txt $direct/$org_type
file copy -force $direct/Scripts/Gen_Pull/Refs/$org_type\_ref.txt $direct/$org_type

cd $direct/$org_type

#Patterns for line searches in the assembly master file ($reflist)
set patb {\n*?.+?}
set patc {.+?\n{1}?}

#Parse the assembly_summary_refseq file for reading.
openfile assembly_summary_refseq.txt
set head_end [string first \n $data]
set reflist [string range $data $head_end end]
#Remove all 'partial' assemblies
regsub -all -line $patb\tPartial$patc $reflist "\n" ass_reflist
set reflist [split [string trim $ass_reflist] \n]


#Open the relevant ref file and parse the species names
openfile $org_type\_ref.txt; 
regsub -all {/} $data {} data
regsub -all {_} $data { } data;
set orglist [split [string trim $data] \n]
regsub -all "{}" $orglist {} orglist
set b 1; #This is a counter for each species being correctly identified in the assembly_summary_refseq
set a 0; #This is a troubleshooting counter

set organism_ids ""
set not_found ""
set multi_error ""

#This section takes each organism name and searches for it in the assembly_summary_refseq.
foreach organ $orglist {
	set match [string trim $organ]

	#This counts the number of hits of the organism name within assembly_summary_refseq, and makes a list from all the options
	set multihits {}
	foreach ref $reflist {
		if {[string first $match $ref] > -1} {
			lappend multihits $ref
		}
	}
	set numbhits [llength $multihits]

	#If the number of hits is one - there is only one assembly version and we use that one regardless of quality.
	if {$numbhits == 1} {
		set ref_line [string trim [lindex $multihits 0]]
		set ids [split $ref_line \t]
		set asm_id [string trim [lindex $ids 0]]
		set asm_name [string trim [lindex $ids 15]]
		append organism_ids "$asm_id\_$asm_name\n"
		incr a
		puts "1 hit"
	#If the no. hits > 1 then we want to find the best quality hit. Primarily this will be either a representative genome (repgen) or a reference genome (refgen)
	} elseif {$numbhits > 1} {
		set x 0
		foreach hit $multihits {
			if {[string first "reference genome" $hit] > -1} {
				set ref_line [string trim $hit]
				set ids [split $ref_line \t]
				set asm_id [string trim [lindex $ids 0]]
				set asm_name [string trim [lindex $ids 15]]
				append organism_ids "$asm_id\_$asm_name\n"
				incr a
				incr x
				puts "Ref hit"
				break
			}
		}
		if {$x == 1} {
		} elseif {$x > 1} {
			puts "ERROR: More than two reference genomes found for $organ"
			error
		} else {
			set x 0
			foreach hit $multihits {
				if {[string first "representative genome" $hit] > -1} {
					set ref_line [string trim $hit]
					set ids [split $ref_line \t]
					set asm_id [string trim [lindex $ids 0]]
					set asm_name [string trim [lindex $ids 15]]
					append organism_ids "$asm_id\_$asm_name\n"
					incr a
					incr x
					puts "Rep hit"
					break
				}
			}
			if {$x == 1} {
			} elseif {$x > 1} {
				puts "ERROR: More than two representative genomes found for $organ"
				error
			#However if there is no refgen, then we want to find the best quality assembly (Complete Genome > Chromosome > Scaffold > Contigs)
			} else {
				set ref_lines {}
				set ref_line ""
				set y 0
				foreach hit $multihits {
					if {[string first "Complete Genome" $hit] > -1 && [string first "representative genome" $hit] == -1} {
						lappend ref_lines $hit
						incr y
					}
				}
				if {[llength $ref_lines] > 0} {
					foreach ref $ref_lines {
						if {[string first "latest" $ref] > 0} {
							set ref_line [string trim $ref]
							set ids [split $ref_line \t]
							set asm_id [string trim [lindex $ids 0]]
							set asm_name [string trim [lindex $ids 15]]
							append organism_ids "$asm_id\_$asm_name\n"
							incr a
							puts "Comp hit"
							break
						}
					}
				} else {
					foreach hit $multihits {
						if {[string first "Chromosome" $hit] > -1 && [string first "representative genome" $hit] == -1} {
							lappend ref_lines $hit
							incr y
						}
					}
					if {[llength $ref_lines] > 0} {
						foreach ref $ref_lines {
							if {[string first "latest" $ref] > 0} {
								set ref_line [string trim $ref]
								set ids [split $ref_line \t]
								set asm_id [string trim [lindex $ids 0]]
								set asm_name [string trim [lindex $ids 15]]
								append organism_ids "$asm_id\_$asm_name\n"
								incr a
								puts "Chr hit"
								break
							}
						}
					} else {
						foreach hit $multihits {
							if {[string first "Scaffold" $hit] > -1 && [string first "representative genome" $hit] == -1} {
								lappend ref_lines $hit
								incr y
							}
						}
						if {[llength $ref_lines] > 0} {
							foreach ref $ref_lines {
								if {[string first "latest" $ref] > 0} {
									set ref_line [string trim $ref]
									set ids [split $ref_line \t]
									set asm_id [string trim [lindex $ids 0]]
									set asm_name [string trim [lindex $ids 15]]
									append organism_ids "$asm_id\_$asm_name\n"
									incr a
									puts "Scaff hit"
									break
								}
							}
						} else {
							foreach hit $multihits {
								if {[string first "latest" $hit] > -1 && [string first "representative genome" $hit] == -1} {
									lappend ref_lines $hit
									incr y
								}
							}
							set ref_line [string trim [lindex $ref_lines 0]]
							set ids [split $ref_line \t]
							set asm_id [string trim [lindex $ids 0]]
							set asm_name [string trim [lindex $ids 15]]
							append organism_ids "$asm_id\_$asm_name\n"
							incr a
							puts "$organ === Contig? hit"
						}
					}
				}
			}
		}
	} else {
		append not_found "$match\n"
	}
	puts "$b/[llength $orglist]"
	incr b
}

#Once we have built a list of the desired organism IDs in the forms GCxxxxx_ASMxxxxx, we output them in the _list.txt file
#The not found organisms (updated assemblies, or other crap) are outputted in the _nf.txt file
regsub -all { } $organism_ids {_} organism_ids
set organism_ids [string trim $organism_ids]
set out [open $direct/$org_type/$org_type\_list.txt a]
set out2 [open $direct/$org_type/$org_type\_nf.txt a]
puts $out $organism_ids
puts $out2 $not_found
close $out
close $out2

###########################################################################
#Here we download the found genomes

#Make a bin folder for the genomes if it doesn't already exist
file mkdir Bin

puts "NAMES EXTRACTED: $a/$b"

#Open the shell file to run the curl and input the correct organism name e.g. Archaea
openfile $direct/Scripts/Gen_Pull/Curl_multi_file.sh
regsub -all {organism} $data $org_type data
regsub -all {filetype} $data $filetype data
set out [open $direct/Scripts/Gen_Pull/Curl_multi_file_$org_type\_$filetype\.sh a]
puts $out $data
close $out

#Give the new .sh file permission to execute and run it
file attributes $direct/Scripts/Gen_Pull/Curl_multi_file_$org_type\_$filetype\.sh -permissions 0777
exec 2>@ stderr $direct/Scripts/Gen_Pull/Curl_multi_file_$org_type\_$filetype\.sh

###########################################################################
#Here we check whether everything was downloaded happily

#Set the downloaded protein names as a list
cd $direct/$org_type/Bin
set protlist [glob *.gz]
cd $direct/$org_type
set i 0
set x 0
set y 0

#Open our _list.txt file full of the desired organism names to download and parse it
openfile $org_type\_list.txt
set namelist [split [string trim $data] \n]
set faillist ""

#If we find that a downloaded file matches one in the _list.txt - then OK
foreach name $namelist {
	if {[regexp $name $protlist] == 1} {
		#puts "ok"
		incr x
	#Else not OK
	} else {
		append faillist "$name\n"
		incr i
	}
	incr y
}

#Create a faillist file to show the organisms that were not downloaded. Then rename the organisms in faillist.txt back to species names.
set faillist [split [string trim $faillist] \n]

set pata {.+?_+?.+?\.+?[0-9]+?}
set patb {.+?\n+?}

set binomial_faillist ""

foreach org $faillist {
	regexp -line $pata $org hit
	regexp -line $hit$patb $ass_reflist match
	set match [split [string trim $match] \t]
	append binomial_faillist "[lindex $match 7]\n"
}

set out [open $direct/$org_type/faillist.txt w]
puts $out [string trim $binomial_faillist]
close $out

###########################################################################

#Get rid of the shell scripts we made to download.
file delete $direct/Scripts/Gen_Pull/Curl_multi_file_$org_type\_$filetype\.sh
#Unzip the Bin directory
cd $direct/$org_type
exec gunzip -r Bin


puts "\nYour downloaded files are in: /Users/(user)/desktop/$org_type"
puts "\nNumber of genome references extracted are: $a out of the refseq list of $b"
puts "\nDownloaded: $x genomes out of a total of $y"
puts "\nThe failed downloads: $i out of $y are stored in /outputdir/faillist.txt"

puts "\nDONE"
return

