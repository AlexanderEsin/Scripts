library(ape)
library(geiger)

setwd("/Users/aesin/Desktop/Mowgli/External_HGT/New_parser")
penalty <- "5"
penalty_refined <- as.vector(read.table(file = paste0("T", penalty, "_refined_groups.tsv"), header = FALSE, sep = "\n")$V1)

tree_tip_number <- vector()
contextual_num = 0
for (group in penalty_refined) {
	tree <- read.tree(file = paste0("/users/aesin/desktop/Mowgli/Final_BS/", group, "_tree.txt"))
	if (length(tree$tip.label) > 49) {
		contextual_num = contextual_num + 1
	}
	tree_tip_number <- c(tree_tip_number, length(tree$tip.label))
}

hist(tree_tip_number, main = paste0("T", penalty, " Refined group sizes"), breaks = seq(from = 1, to = max(tree_tip_number)+50, by = 50))