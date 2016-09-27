library(dplyr)
library(ggplot2)

setwd("/users/aesin/desktop/Geo_analysis/Geo_omes")

x <- read.delim("Geobac_num_proteins.tsv")
x_sorted <- arrange(x, Number.of.genes)
x_sorted$Name <- reorder(x_sorted$Name, x_sorted$Number.of.genes)

p1 <- ggplot(x_sorted, aes(x = Name, y = Number.of.genes)) + geom_point(stat = "identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1, colour = "black", size = 12), axis.title.x = element_blank(), axis.text.y = element_text(colour = "black")) + ylab("Number of predicted proteins") + ggtitle("Number of predicted proteins across Geobacillus taxa") + ylim (2000, 4000)

ggsave(p1, file="Geobac_protein_distrib.png")