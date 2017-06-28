### General util procs ###

############################################################
## General file manipulation ##
# Open and read in a file #
proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return $data
}

proc multiputs {args} {
    if { [llength $args] == 0 } {
        error "Usage: multiputs ?channel ...? string"
    } elseif { [llength $args] == 1 } {
        set channels stdout
    } else {
        set channels [lrange $args 0 end-1]
    }
    set str [lindex $args end]
    foreach ch $channels {
        puts $ch $str
    }
}

############################################################
## String formatting ##
# The var is the inserted number variable. {num 3} gives the option for separating larger intervals with a default of 3. {char ,} gives the option of a different separating character, with the default being a comma. #
proc commas {var {num 3} {char ,}} {
    set len [string length $var]
    set first [expr $len - $num]
    set x {}
    while {$len > 0} {
        # grab left num chars
        set lef [string range $var $first end] 
        if {[string length $x] > 0} {
            set x "${lef}$char${x}"
        } else {
            set x $lef
        }
        # grab everything except left num chars
        set var [string range  $var 0 [expr $first -1]]
        set len [string length $var]
        set first [expr {$len - $num}]
    }
    return $x
}

# Escape all the special characters in the variable and return the escaped string #
proc escape_special {var} {
  set _var [string map [list \[ {\[} \] {\]} \* {\*} \? {\?} \\ {\\}] $var]
  return $_var
}

# Take a string and reverse the order of each character #
proc string_reverse {some_string} {
  set string_split [split $some_string {}]
  set reversed_list [lreverse $string_split]
  
  set reversed_string [join $reversed_list {}]
  return $reversed_string
}

############################################################
## Fasta specific ##

# RETIRED: 
# This proc takes the read input from a typical fasta files and splits it into a list of genes #
# proc split_genes {fasta} {
#   global genes
#   regsub -all {>} [string trim $fasta] {£>} fasta
#   set fasta [split $fasta £]
#   regsub -all "{}" $fasta {} genes
#   return $genes
# }

# This does the same as above but x15 faster #
proc split_genes {fasta} {
  global genes

  set temp [split [string range [string trim $fasta] 1 end] >]
  set genes {}
  foreach seq $temp {
    lappend genes ">$seq"
  }
  return $genes
}

# Format raw sequence (no header) into a fasta-like linebroken style #
proc linebreak {s {width 80}} {
   global res
   set res ""
   while {[string length $s]>$width} {
       set res "$res[string range $s 0 [expr $width - 1]]\n"
       set s [string range $s $width end]
   }
   set res "$res$s"
   return $res
}

# As above but with a width of 60 - used in this format in Geo_v_all/11.Finalise_cores_split_for_trees #
proc linebreak_align {s {width 60}} {
   global res
   set res ""
   while {[string length $s]>$width} {
       set res "$res[string range $s 0 59]\n"
       set s [string range $s 60 end]
   }
   set res "$res$s"
   return $res
}

############################################################
## List manipulation ##
# Order a list in reverse dictionary order #
proc reverse_dict {lst} {
    global revdict
    set revdict {}
    set lst [lsort -dictionary $lst]
    while {[llength $lst] > 0} {
        set element [lindex $lst end]
        lappend revdict $element
        set lst [lrange $lst 0 end-1]
    }
    return $revdict
}

# Find and return duplicates in any given list #
proc dups xs {
    set dups {}
    set xs [lsort $xs]
    set len [expr {[llength $xs]-1}]
    for {set i 0} {$i < $len} {} {
        set x [lindex $xs $i]
        if {$x == [lindex $xs [incr i]]} {
            lappend dups $x
            # discard multiple repeats
            while {[lindex $xs $i] == $x} { incr i }
        }
    }
    return $dups
}

# Sum all the elements in a list #
proc ladd {l} {::tcl::mathop::+ {*}$l}

# Remove every copy of element from list #
proc lremove {list element} {
    set idx_list [lsearch -all $list $element]
    foreach idx [lreverse $idx_list] {
        set list [lreplace $list $idx $idx]
    }
    return $list
}

# Remove any empty elements in a list #
proc lremove_empty {list} {
  set idx_list [lsearch -all $list {}]
  foreach idx [lreverse $idx_list] {
      set list [lreplace $list $idx $idx]
  }
  return $list
}

# Count elements in a list #
proc lcount list {
  foreach x $list {lappend arr($x) {}}
  set res {}
  foreach name [array names arr] {
    lappend res [list $name [llength $arr($name)]]
  }
  return $res
}

############################################################
## Maths ##
# Rounding proc to two dp #
proc tcl::mathfunc::roundto {value {decimalplaces 2}} {
  expr {round(10**$decimalplaces*$value)/10.0**$decimalplaces}
}

proc RandomInt {min max} {
    return [expr {int(rand()*($max-$min+1)+$min)}]
}

############################################################
## Progress meters ##
proc spinner {} {
    global spinnerIdx
    if {[incr spinnerIdx] > 3} {
        set spinnerIdx 0
    }
    set spinnerChars {/ - \\ |}
    puts -nonewline "\b[lindex $spinnerChars $spinnerIdx]"
    flush stdout
}

# Progress BAR #
proc progress_init {tot} {
   set ::progress_start     [clock seconds]
   set ::progress_last      0
   set ::progress_last_time 0
   set ::progress_tot       $tot
}

# We update if there's a 5% difference or a 5 second difference

proc progress_tick {cur} {
   set now [clock seconds]
   set tot $::progress_tot

   if {$cur > $tot} {
       set cur $tot
   }
   if {($cur >= $tot && $::progress_last < $cur) ||
       ($cur - $::progress_last) >= (0.05 * $tot) ||
       ($now - $::progress_last_time) >= 5} {
       set ::progress_last $cur
       set ::progress_last_time $now
       set percentage [expr round($cur*100/$tot)]
       set ticks [expr $percentage/2]
       if {$cur == 0} {
           set eta   ETA:[format %7s Unknown]
       } elseif {$cur >= $tot} {
           set eta   TOT:[format %7d [expr int($now - $::progress_start)]]s
       } else {
           set eta   ETA:[format %7d [expr int(($tot - $cur) * ($now - $::progress_start)/$cur)]]s
       }
       set lticks [expr 50 - $ticks]
       set str "[format %3d $percentage]%|[string repeat = $ticks]"
       append str "[string repeat . $lticks]|[format %8d $cur]|$eta\r"
       puts -nonewline stdout $str
       if {$cur >= $tot} {
           puts ""
       }
       flush stdout
   }
}

############################################################
## Compare two files ##
proc cmpFilesChunked {file1 file2 {chunksize 16384}} {
  set f1 [open $file1]; fconfigure $f1 -translation binary
  set f2 [open $file2]; fconfigure $f2 -translation binary
  while {1} {
    set d1 [read $f1 $chunksize]
    set d2 [read $f2 $chunksize]
    set diff [string compare $d1 $d2]
    if {$diff != 0 || [eof $f1] || [eof $f2]} {
      close $f1; close $f2
      return $diff
    }
  }
}

############################################################