source ~/Documents/Scripts/General_utils.tcl

cd ~/Documents/Scripts

proc findFiles {basedir extension} {

    # Fix the directory name, this ensures the directory name is in the
    # native format for the platform and contains a final directory seperator
    set basedir [string trimright [file join [file normalize $basedir] { }]]
    set fileList {}

    # Look in the current directory for matching files, -type {f r}
    # means ony readable normal files are looked at, -nocomplain stops
    # an error being thrown if the returned list is empty
    foreach fileName [glob -nocomplain -type f -path $basedir *$extension] {
        lappend fileList $fileName
    }

    # Now look for any sub direcories in the current directory
    foreach dirName [glob -nocomplain -type {d  r} -path $basedir *] {
        # Recusively call the routine on the sub directory and append any
        # new files to the results
        set subDirList [findFiles $dirName $extension]
        if { [llength $subDirList] > 0 } {
            foreach subDirFile $subDirList {
                lappend fileList $subDirFile
            }
        }
    }
    return $fileList
 }

set all_tcl_scripts [findFiles /Users/aesin/Documents/Scripts/ .tcl]

foreach script $all_tcl_scripts {
	set script_data [openfile $script]
	if {[regexp -all {Inside_group_tips} $script_data] > 0} {
		puts $script
	}
}