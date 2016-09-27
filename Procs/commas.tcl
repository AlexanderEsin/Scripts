# The var is the inserted number variable. {num 3} gives the option for separating larger intervals with a default of 3. {char ,} gives the option of a different separating character, with the default being a comma. #
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