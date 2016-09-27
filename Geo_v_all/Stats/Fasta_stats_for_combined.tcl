###########################################################################
### PROCS ###
source ~/Dropbox/Scripts/General_utils.tcl
###########################################################################


set ival 2.0
set direct /users/aesin/desktop/Geo_analysis/Geo_v_all/$ival

###########################################################################

cd $direct
file mkdir Geo_only_groups/Group_fastas
file mkdir Anoxy_geo_only_groups/Group_fastas

cd $direct/Final_family_groups

set fasta_files [glob *.faa]
set fasta_files [lsort -dictionary $fasta_files]

set output_table {}
lappend output_table "Group\tTotal proteins\tTotal Geobacillus"
set 19_geo {}
set 19_geo_only {}
set 19_geo_more_than_context {}

set i 0
set 19_geo_counter 0
set geo_only_counter 0
set geo_anoxy_only 0
set 19_geo_more_than_context_counter 0
set total_prot_master 0
set total_geobac_master 0


foreach file $fasta_files {
    openfile $file
    set group_no [string range $file 0 end-4]

    set total_prot_number [regexp -all {>} $data]

    set total_prot_master [expr $total_prot_master + $total_prot_number]

    if {$total_prot_number < 4} {
        continue
    }

    set total_geobac_number [regexp -all {\(Geobac\)} $data]

    set total_anoxy_number [regexp -all {Anoxybacillus} $data]

    set total_geobac_master [expr $total_geobac_master + $total_geobac_number]

    if {$total_geobac_number == 19} {
        incr 19_geo_counter
        lappend 19_geo $group_no
    }

    set total_anoxy_geo_number [expr $total_geobac_number + $total_anoxy_number]

    if {$total_geobac_number == $total_prot_number} {
        incr geo_only_counter
        lappend geo_only "$group_no\t$total_geobac_number"
        file copy -force $file $direct/Geo_only_groups/Group_fastas
    } elseif {$total_anoxy_geo_number == $total_prot_number} {
        incr $geo_anoxy_only
        file copy -force $file $direct/Anoxy_geo_only_groups/Group_fastas
    }

    if {$total_geobac_number == 19 && [expr $total_geobac_number * 2] > $total_prot_number} {
        incr 19_geo_more_than_context_counter
        lappend 19_geo_more_than_context $group_no
    }

    lappend output_table "$group_no\t$total_prot_number\t$total_geobac_number"
    incr i
    puts "Counting.... $i/[llength $fasta_files]"
}

set geo_only_string [join $geo_only \n]

set out [open $direct/Stats/Combined_group_fastas_stats.tsv w]
puts $out [join $output_table \n]
close $out

set output_txt "No. groups total: $i\nNo. proteins total: $total_prot_master\nNo. Geobacilli total: $total_geobac_master\nNo. groups with exactly 19 Geobacilli: $19_geo_counter\nNo. groups with only Geobacilli: $geo_only_counter\nNo. groups with exactly 19 Geobacilli and fewer than 19 other taxa: $19_geo_more_than_context_counter"
set out [open $direct/Stats/Combined_group_fastas_log.txt w]
puts $out $output_txt
close $out

puts $output_txt





