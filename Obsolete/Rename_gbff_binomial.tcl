#!/usr/local/bin/tclsh

## Should translate the ncbi ids used by Refseq into the binomial names - to match those used by the proteome files ##
## Output is a replacement of the ncbi_id named file with the binomial named file ##
## Additionally, a translation table is made for easier reference downstream to link the ncbi id and the binomial name ##

## Table looks like this:
# GCF_000003925.1_ASM392v1	Bacillus_mycoides_DSM_2048
# GCF_000005825.2_ASM582v2	Bacillus_pseudofirmus_OF4
# GCF_000006175.1_ASM617v2	Methanococcus_voltae_A3
# GCF_000006605.1_ASM660v1	Corynebacterium_jeikeium_K411
# GCF_000006625.1_ASM662v1	Ureaplasma_parvum_serovar_3_str._ATCC_700970

###########################################################################

## Source procs ##
source ~/Dropbox/Scripts/General_utils.tcl

proc taxa_name_clean {name} {
	global newname
	regsub -all "/" $name "_" newname
	regsub -all -- {\-} $newname "_" newname
	regsub -all {=} $newname {} newname
	regsub -all "  " $newname " " newname
	regsub -all " " $newname "_" newname
	return $newname
}

###########################################################################

## Set variables ##
set direct /users/aesin/desktop/Clean_Genomes
set folder Prok_gbff_files

set file_ext "_genomic.gbff"
set output_table_switch 1

## Open the assembly reference file ##

cd ~/Dropbox/Scripts/Gen_Pull
openfile assembly_summary_refseq.txt
set reflist $data

## Get a list of the gbff files ##

cd $direct/$folder
puts "Processing files from [exec pwd]"
set gbff_files [glob *$file_ext]

## Process and rename each file ##

set i 1

foreach gbff_f $gbff_files {
	regsub "$file_ext" $gbff_f {} old_name

	openfile $gbff_f
	set org_name_line [lindex [regexp -line -inline {ORGANISM.+} $data] 0]
	set org_name [string trim [string range $org_name_line [string first " " $org_name_line] end]]

	set new_name [taxa_name_clean $org_name]
	file rename -force $gbff_f $new_name$file_ext

	lappend translation_table "$old_name\t$new_name"

	puts "Renamed ... $i/[llength $gbff_files]"
	incr i
}

if {$output_table_switch == 1} {
	set out [open $direct/$folder\_ID_binomial_translation.tsv w]
	puts $out [join $translation_table \n]
	close $out
}