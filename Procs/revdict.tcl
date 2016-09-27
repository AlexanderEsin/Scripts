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
