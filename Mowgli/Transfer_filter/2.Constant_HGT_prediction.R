#!/usr/local/bin/Rscript

## Compare the results values at different transfer penalty costs and extract those groups which are positive for HGT signal into Anoxy/Geo at all the penalties ##
library(stringr)
library(plyr)
library(dplyr)
library(reshape2)
library(lazyeval)
library(gtools)

ReadInRawTrans	<- function(penalty) {
	# Find the file
	file_name	<- paste0("T", penalty, "_D2_L1_results.txt")

	# Transfer data
	trans_data	<- read.table(file = paste0(in_path, file_name), header = FALSE, sep = "\t")
	colnames(trans_data)	<- c("Group", paste0("HGT_", penalty), paste0("AGRoot_", penalty))
	return(trans_data)
}

penalty_l <- list(3, 4, 5, 6)

#######
## Get and order the results files by dictionary order (package = gtools) ##
in_path		<- "/Users/aesin/Desktop/Mowgli/Raw_predictions/"
out_path	<- "/Users/aesin/Desktop/Mowgli/Constant_events/"
dir.create(out_path, recursive = TRUE, showWarnings = FALSE)


# Get raw transfer prediction results files
raw_results_data_l	<- lapply(penalty_l, ReadInRawTrans)
name_vector			<- unlist(lapply(penalty_l, function(penalty) paste0("T", penalty)))
names(raw_results_data_l)	<- name_vector

temp	<- lapply(penalty_l, function(transfer_penalty) {

	hgt_pen_l	<- penalty_l[penalty_l <= transfer_penalty]
	ver_pen_l	<- penalty_l[penalty_l >= transfer_penalty]

	hgt_pen_colnames	<- unlist(lapply(hgt_pen_l, function(penalty) paste0("T", penalty)))
	ver_pen_colnames	<- unlist(lapply(ver_pen_l, function(penalty) paste0("T", penalty)))

	hgt_data	<- raw_results_data_l[hgt_pen_colnames]
	ver_data	<- raw_results_data_l[ver_pen_colnames]

	hgt_merged	<- Reduce(function(df1, df2) inner_join(df1, df2, by = "Group"), hgt_data)
	ver_merged	<- Reduce(function(df1, df2) inner_join(df1, df2, by = "Group"), ver_data)

	hor_conditions	<- paste(colnames(hgt_merged)[c(seq(2, ncol(hgt_merged), by = 2))], "== 1", collapse = " & ")
	ver_conditions	<- paste(colnames(ver_merged)[c(seq(2, ncol(ver_merged), by = 2))], "== 0", collapse = " & ")

	# Root filter
	hor_root_filter	<- paste(colnames(hgt_merged)[c(seq(3,ncol(hgt_merged), by = 2))], "== FALSE", collapse = " & ")
	ver_root_filter	<- paste(colnames(ver_merged)[c(seq(3,ncol(ver_merged), by = 2))], "== FALSE", collapse = " & ")

	true_hor_intersect	<- filter_(filter_(hgt_merged, hor_conditions), hor_root_filter)
	true_ver_intersect	<- filter_(filter_(ver_merged, ver_conditions), ver_root_filter)


	hor_pen_outname_str	<- paste0("HGT_const_t", paste(unlist(hgt_pen_l), collapse = "_t"), ".tsv")
	ver_pen_outname_str	<- paste0("Ver_const_t", paste(unlist(ver_pen_l), collapse = "_t"), ".tsv")


	write.table(true_hor_intersect$Group, file = paste0(out_path, hor_pen_outname_str), sep = "\n", quote = FALSE, col.names = FALSE, row.names = FALSE)
	write.table(true_ver_intersect$Group, file = paste0(out_path, ver_pen_outname_str), sep = "\n", quote = FALSE, col.names = FALSE, row.names = FALSE)

})