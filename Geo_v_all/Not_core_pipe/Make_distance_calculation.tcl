set direct /users/aesin/Desktop/Geo_v_all/-100
set rbbh_dir Rbbh_all_-100
set hitfile RBBH_Hitfile_-100.txt

proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}

proc commas {var {num 3} {char ,}} {
    set len   [string length $var]
    set first [expr $len - $num]
    set x     {}
    while {$len > 0} {
        # grab left num chars
        set lef [string range $var $first end] 
        if {[string length $x] > 0} {
            set x   "${lef}$char${x}"
        } else {
            set x   $lef
        }
        # grab everything except left num chars
        set var [string range  $var 0 [expr $first -1]]
        set len   [string length $var]
        set first [expr {$len - $num}]
    }
    return $x
}


package require sqlite3

###########################################################################

namespace import ::tcl::mathop::*
proc average list {expr {[+ {*}$list]/double([llength $list])}}

cd /users/aesin/desktop/Clean_Proteomes
sqlite3 db1 all_prot_geo_db

cd $direct

openfile Group_list_2.5.tsv

set groups [split [string trim $data] \n]
set groups [lrange $groups 141 141]

set master_table "Group_no\tGeobac_dist"

set header_counter 0
set i 0
set z 1
set outliers {}
set master_list ""
foreach group $groups {
    # There have two be at least two Geobacilli in the group for it to be evaluated #
    if {[regexp -all {\(Geobac\)} $group] < 2} {
        incr z
        continue
    }
    set genes [split [string trim $group] \t]
    set group_no [string trim [lindex $genes 0]]
    set genes [lrange $genes 1 end]

    if {$header_counter == 0} {
        set a 1
        foreach gene $genes {
            set master_table "$master_table\tGene_$a"
            incr a
        }
        incr header_counter
    }

    set geobacs {}
    # Get a list of all the geobacilli ids in the group #
    foreach gene $genes {
        if {[string first "(Geobac)" $gene] > 0} {
            set id [string trim [string range $gene 0 [string first " " $gene]]]
            lappend geobacs $id
        }
    }
    # If the group only contains Geobacilli, the group is ignored #
    set geobac_no [llength $geobacs]
    if {$geobac_no == [llength $genes]} {
        incr z
        continue
    }

    ### Get the average distance between the Geobacilli (based on eval) ###
    set geobac_eval_list {}
    set geobac_total_eval 0

    foreach geobac_x $geobacs {
        # For each Geobacillus, extract the file_name #
        db1 eval {SELECT file_name FROM t1 WHERE id = $geobac_x} {
            set name_x "$file_name"
        }
        set name_x [string range $name_x 0 end-4]
        foreach geobac_y $geobacs {
            # Do not compare the id to itself #
            if {$geobac_x != $geobac_y} {
                # Extract the file_name for the second Geobacillus #
                db1 eval {SELECT file_name FROM t1 WHERE id = $geobac_y} {
                    set name_y "$file_name"
                }
                set name_y [string range $name_y 0 end-4]
                # If the two geobac ids are from the same organism (potential paralogs), ignore that interaction since we do not have a self_v_self blast #
                if {$name_x != $name_y} {
                    # Pick the correct file name, since only one of RBBH_x&y.txt or RBBH_y&x.txt will exist #
                    set fname "RBBH_$name_x\&$name_y\.txt"
                    if {[file exists $direct/$rbbh_dir/Bin/$fname] == 0} {
                        set fname "RBBH_$name_y\&$name_x\.txt"
                    }

                    openfile $direct/$rbbh_dir/Bin/$fname
                    set data "\n$data\n"
                    # Use string first rather than regular expressions to find the correct entry and get the evalue. There is a possibility that two Geobacilli may not have hit each other at the blast level but were grouped together later based on similar hits to other taxa #
                    if {[string first "\n$geobac_x" $data] > 0} {
                        set eval "ERROR"
                        set temp_string "[string trim [string range $data [string first "\n$geobac_x" $data] [string first "\n$geobac_x" $data]+50]]\n"
                        set hit_string [string trim [string range $temp_string 0 [string first \n $temp_string]]]
                        set eval [lindex [split $hit_string \t] 2]
                        if {$eval == 0.0} {
                            set eval 1e-200
                        }
                        set eval [expr log10($eval)]
                        puts "$geobac_x === $geobac_y === $eval"
                        lappend geobac_eval_list $eval

                        set geobac_total_eval [expr $geobac_total_eval + $eval]

                        #puts "$fname === $geobac_x === $geobac_y === $eval"

                    }                   
                }
            }
        }
    }
    if {$geobac_total_eval == 0} {
        incr z
        continue
    }
    set deltaGeo [expr [average $geobac_eval_list] + 200]
    set group_table "Group_$group_no\t$deltaGeo"
    
    ## For each gene get the average against every geobac ##
    set y 1
    foreach gene $genes {
        set gene [string trim [string range $gene 0 [string first " " $gene]]]
        set gene_eval_list {}
        set gene_total_eval 0
        # For each id, extract the file_name #
        db1 eval {SELECT file_name FROM t1 WHERE id = $gene} {
            set name_x "$file_name"
        }
        set name_x [string range $name_x 0 end-4]
        if {[string first "Geobacillus" $name_x] > -1} {
            puts "IGNORE ========= $name_x"
            incr aa
            continue
        } else {
            foreach geobac $geobacs {
                # Do not compare the id to itself #
                if {$gene != $geobac} {
                    # Extract the file_name for the second Geobacillus #
                    db1 eval {SELECT file_name FROM t1 WHERE id = $geobac} {
                        set name_y "$file_name"
                    }
                    set name_y [string range $name_y 0 end-4]
                    # If the two geobac ids are from the same organism (potential paralogs), ignore that interaction since we do not have a self_v_self blast #
                    if {$name_x != $name_y} {
                        # Pick the correct file name, since only one of RBBH_x&y.txt or RBBH_y&x.txt will exist #
                        set fname "RBBH_$name_x\&$name_y\.txt"
                        if {[file exists $direct/$rbbh_dir/Bin/$fname] == 0} {
                            set fname "RBBH_$name_y\&$name_x\.txt"
                        }

                        openfile $direct/$rbbh_dir/Bin/$fname
                        set data "\n$data\n"
                        # Use string first rather than regular expressions to find the correct entry and get the evalue. There is a possibility that a gene may not have hit a specific Geobacillus at the blast level but were grouped together later based on similar hits to other taxa #
                        if {[string first "\n$gene" $data] > 0} {
                            set eval "ERROR"
                            set temp_string "[string trim [string range $data [string first "\n$gene" $data] [string first "\n$gene" $data]+50]]\n"
                            set hit_string [string trim [string range $temp_string 0 [string first \n $temp_string]]]
                            set eval [lindex [split $hit_string \t] 2]
                            if {$eval == 0.0} {
                                set eval 1e-200
                            }
                            set eval [expr log10($eval)]
                            # puts "$geobac === $eval"
                            lappend gene_eval_list $eval

                            set gene_total_eval [expr $gene_total_eval + $eval]
                        }                   
                    }
                }
            }
            if {$gene_total_eval != 0} {
                set deltaGene [expr [average $gene_eval_list] + 200]
                if {$deltaGene < $deltaGeo} {
                    lappend outliers $gene                }
                set group_table "$group_table\t$deltaGene"
            }
            puts "$z / [llength $groups] ==== $y / [llength $genes]"
            incr y
        }
    }
    incr z
    set master_table "$master_table\n$group_table"

    foreach outlier $outliers {
        set hit_counter 0
        set gene $outlier
        set gene_eval_list {}
        set gene_total_eval 0
        # For each id, extract the file_name #
        db1 eval {SELECT file_name FROM t1 WHERE id = $gene} {
            set name_x "$file_name"
        }
        db1 eval {SELECT binomial FROM t1 WHERE id = $gene} {
            set binomial_outlier "$binomial"
        }
        set binomial [string range $binomial 1 end-1]
        set name_x [string range $name_x 0 end-4]
        if {[string first "Geobacillus" $name_x] > -1} {
            puts "IGNORE ========= $name_x"
            incr aa
            continue
        } else {
            foreach geobac $geobacs {
                # Do not compare the id to itself #
                if {$gene != $geobac} {
                    # Extract the file_name for the second Geobacillus #
                    db1 eval {SELECT file_name FROM t1 WHERE id = $geobac} {
                        set name_y "$file_name"
                    }
                    db1 eval {SELECT binomial FROM t1 WHERE id = $geobac} {
                        set binomial_geo "$binomial"
                    }
                    set name_y [string range $name_y 0 end-4]
                    # If the two geobac ids are from the same organism (potential paralogs), ignore that interaction since we do not have a self_v_self blast #
                    if {$name_x != $name_y} {
                        # Pick the correct file name, since only one of RBBH_x&y.txt or RBBH_y&x.txt will exist #
                        set fname "RBBH_$name_x\&$name_y\.txt"
                        if {[file exists $direct/$rbbh_dir/Bin/$fname] == 0} {
                            set fname "RBBH_$name_y\&$name_x\.txt"
                        }

                        openfile $direct/$rbbh_dir/Bin/$fname
                        set data "\n$data\n"
                        # Use string first rather than regular expressions to find the correct entry and get the evalue. There is a possibility that a gene may not have hit a specific Geobacillus at the blast level but were grouped together later based on similar hits to other taxa #
                        if {[string first "\n$gene" $data] > 0} {
                            set eval "ERROR"
                            set temp_string "[string trim [string range $data [string first "\n$gene" $data] [string first "\n$gene" $data]+50]]\n"
                            set hit_string [string trim [string range $temp_string 0 [string first \n $temp_string]]]
                            set eval [lindex [split $hit_string \t] 2]
                            if {$eval == 0.0} {
                                set eval 1e-200
                            }
                            set eval [expr log10($eval)]
                            puts "$geobac === $binomial_geo === $eval"
                            incr hit_counter
                            lappend gene_eval_list $eval

                            set gene_total_eval [expr $gene_total_eval + $eval]
                        } else {
                            puts "MISSING == $geobac === $binomial_geo"
                        }                   
                    }
                }
            }
        }
        set deltaGene [expr [average $gene_eval_list] + 200]
        puts "\n\n$outlier == $binomial_outlier == $deltaGene === $hit_counter\n\n"
    }
}

set out [open $direct/Master_delta_geo_142.tsv w]
puts $out [string trim $master_table]
close $out

puts [llength $geobacs]
puts $outliers

