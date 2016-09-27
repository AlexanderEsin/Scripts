proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}
