#!/usr/bin/env Rscript

# Any libraries explicitly used in the script
require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load(tidyverse, fs)

group		<- "Bacillus"
master_path	<- file.path("/home/ade110/Work", group)
mowOut_path	<- file.path(master_path, "Mowgli", "Mowgli_output")

penalty_dirs	<- list.dirs(mowOut_path, full.names = TRUE, recursive = FALSE)

result_list	<- lapply(penalty_dirs, function(pen_dir) {
	message(pen_dir)

	pen_value	<- str_extract(pen_dir, "(?<=_)[0-9]+")
	subdirs		<- list.dirs(pen_dir, full.names = TRUE, recursive = FALSE)
	message(subdirs[1:2])
	completed_l	<- lapply(subdirs, function(subdir) file.exists(file.path(subdir, "reconciled.ok")))

	numCompleted	<- length(which(unlist(completed_l)))
	outTbl	<- tibble(Penalty = pen_value, numCompleted = numCompleted)
	return(outTbl)
})

result_df	<- bind_rows(result_list)
write_tsv(result_df, path = file.path(mowOut_path, "completedRecon.tsv"))