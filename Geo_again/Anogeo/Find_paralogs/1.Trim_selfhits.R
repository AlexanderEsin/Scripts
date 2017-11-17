#!/usr/local/bin/Rscript

## First argument is file path to Blast output file, second argument is the output directory
args <- commandArgs(trailingOnly = TRUE)

input_path	<- as.character(args[1])
output_path	<- as.character(args[2])
file_name	<- basename(input_path)

## Read in the Blast output file
para_blast_tbl	<- read.table(input_path, header = FALSE, sep = "\t")

## Remove against self hits
para_blast_trim	<- para_blast_tbl[-which(para_blast_tbl$V1 == para_blast_tbl$V2),]

## For the remainder, pick only the top blast hit
para_blast_top	<- para_blast_trim[!duplicated(para_blast_trim$V1),]

write.table(file = paste0(output_path, "/", file_name), x = para_blast_top, sep = "\t", quote = FALSE, row.name = FALSE, col.names = FALSE)