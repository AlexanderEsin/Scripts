## Use this script to find discrepancies in the number of Geobacilli between the final (combined) group families and those found at the lowest stringnecy evalue (-10), output those groups that demonstrate a different geobac number into a table with the group number and the number of missing geobacs. See Group_size_distrib.R for a plot script for that table ##

source ~/Dropbox/Scripts/General_utils.tcl

set ival 2.0
set direct /users/aesin/desktop/Geo_analysis/Geo_v_all/$ival

###########################################################################

set combined_groups $direct/Final_family_groups
set 10_eval_groups $direct/-10/Group_fastas_$ival

cd $combined_groups
set fasta_files [glob *.faa]
set fasta_files [lsort -dictionary $fasta_files]

set output_table {}
set header "Group number\tNumber of missing Geobacilli sequences"
lappend output_table $header

set i 0
set geo_missing_counter 0
set total_comb_geo_number 0
set total_10_geo_number 0

foreach file $fasta_files {
    openfile $file
    set combined_genes $data
    set group_no [string range $file 0 end-4]

    # Get the number of geobac genes in the Final_family_groups fasta file (aka combined group) #
    set combined_geobac_number [regexp -all {\(Geobac\)} $combined_genes]
    set total_comb_geo_number [expr $total_comb_geo_number + $combined_geobac_number]

    # Get the number of geobac genes in the lowest eval stringnecy fasta file (aka -10 eval fasta file) #
    openfile $10_eval_groups/$file
    set 10_eval_genes $data
    set 10_eval_geobac_number [regexp -all {\(Geobac\)} $10_eval_genes]
    set total_10_geo_number [expr $total_10_geo_number + $10_eval_geobac_number]

    # If the two geobac numbers are not the same ... #
    if {$combined_geobac_number != $10_eval_geobac_number} {
        set no_missing_geobacs [expr $10_eval_geobac_number - $combined_geobac_number]
        set entry "$group_no\t$no_missing_geobacs"
        lappend output_table $entry
        incr geo_missing_counter
    }

    incr i
    puts "Counting.... $i/[llength $fasta_files]"
}

set output_table_str [join $output_table \n]
set out [open $direct/Stats/Missing_geobac_families.tsv w]
puts $out $output_table_str
close $out

puts "\n$output_table_str"
puts "\n=====Numbers below should match the stats in Geo_v_all/2.0/Stats/Combined_v_-10_groups.txt====="
puts "Total number of Geobacilli in combined set: $total_comb_geo_number"
puts "Total number of Geobacilli in the -10 eval set: $total_10_geo_number"





