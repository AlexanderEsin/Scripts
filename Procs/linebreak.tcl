#This proc remakes the 80-wide fasta entries from a flat CDS sequence

proc linebreak {s {width 80}} {
   global res
   set res ""
   while {[string length $s]>$width} {
       set res "$res[string range $s 0 [expr $width - 1]]\n"
       set s [string range $s $width end]
   }
   set res "$res$s"
}