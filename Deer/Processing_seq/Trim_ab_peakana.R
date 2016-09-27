library(stringr)

setwd("/Volumes/molsystems$/Deer_project/Sequences/6113_seq_13.8.15/Ab1_peak_reports_raw")

for (ab_file in dir()) {
  x<-read.csv(file = ab_file, skip = 15, header=T, stringsAsFactors=F)
  x_new<-x[-which(x$BaseCall=="-"),]
  x_new_2<-x_new[which(x_new$QualityValue >= 10),]
  x_new_3<-x_new_2[which(nchar(x_new_2$SequenceWindow) == 7),]
  x_new_4<-x_new_3[complete.cases(x_new_3$MaxSig7Scan_Filtered_G_Ratio),]
  trim_suff<-str_sub(ab_file, 0, -5)
  output_name<-paste("Stripped_", trim_suff, ".tsv", sep = "")
  write.table(x_new_4, file = output_name, sep = "\t", col.names=T)
}

