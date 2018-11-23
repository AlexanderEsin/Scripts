library(stringr)
library(dplyr)
library(ggplot2)


setwd("~/Desktop/doric6.5/")
oricBactDat <-
	read.table(
		"bacteria_record.dat",
		sep = "\t",
		header = TRUE,
		stringsAsFactors = FALSE
	)

oriC_data 	<-
	str_split(
		str_replace(
			str_trim(oricBactDat$oriC.location),
			pattern = " nt[\\*]?",
			replacement = ""
		),
		pattern = "\\.\\.",
		simplify = T
	)

oriC_df		<- data.frame(
	dORI_AC = oricBactDat$DoriC.AC,
	organism = oricBactDat$Organsim,
	genomeLength = as.numeric(str_trim(
		str_replace(
			oricBactDat$Genome.Length,
			pattern = " nt[\\*]?",
			replacement = ""
		)
	)),
	oriC_start = as.numeric(oriC_data[, 1]),
	oriC_end = as.numeric(oriC_data[, 2]),
	stringsAsFactors = FALSE
)


dnaA_multifilter	<-
	subset(oricBactDat[!str_detect(oricBactDat$dnaA.location, ",|-"), ], select = c(DoriC.AC, dnaA.location))
dnaA_data 			<-
	str_split(
		str_replace(
			str_trim(dnaA_multifilter$dnaA.location),
			pattern = " nt[\\*]?",
			replacement = ""
		),
		pattern = "\\.\\.",
		simplify = T
	)

dnaA_df				<- data.frame(
	dORI_AC = dnaA_multifilter$DoriC.AC,
	dnaA_start = as.numeric(dnaA_data[, 1]),
	dnaA_end = as.numeric(dnaA_data[, 2]),
	stringsAsFactors = FALSE
)


combined_df			<- left_join(dnaA_df, oriC_df, by = "dORI_AC")
combined_df$dnaRel	<-	combined_df$dnaA_start / combined_df$genomeLength
combined_df$oriRel	<-	combined_df$oriC_start / combined_df$genomeLength
ori_dnaA_dist		<-	signif(abs(combined_df$dnaRel - combined_df$oriRel), digits = 3)
combined_df$dist	<-	ifelse(ori_dnaA_dist > 0.5, 1 - ori_dnaA_dist, ori_dnaA_dist)

AG_positions		<-	combined_df[grep(x = combined_df$organism, pattern = "Geobacillus|Anoxybacillus"), ]


ggplot(data = combined_df, aes(x = dist, y = ..scaled..)) + stat_density(geom = "line") + theme_classic()
