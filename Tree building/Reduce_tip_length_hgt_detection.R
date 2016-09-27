### Reducing the lengths of tips - for the hgt bipartitions method. ###

# perl /Users/aesin/Downloads/hgt-detection_3.22/run_hgt.pl -inputfile Input.txt


library(stringr)
library(ape)

setwd("/users/aesin/desktop/Bipartition_test")

reference_tree <- "Test_adding_branch_len.txt"
reference_data <- read.tree(reference_tree)

for (n in 1:length(reference_data$tip.label)) {

	tip_name <- reference_data$tip.label[n]
	tip_name_len <- nchar(tip_name)
	if (tip_name_len > 48) {
		new_tip_name <- str_sub(tip_name, 1, 48)
		reference_data$tip.label[n] <- new_tip_name
	}
}

write.tree(reference_data, file = "Reference_tree_trim_long_tips.txt")


setwd("/users/aesin/desktop/Bipartition_test")

reference_tree <- "Relabelled_656_nbs.txt"
reference_data <- read.tree(reference_tree)

for (n in 1:length(reference_data$tip.label)) {

	tip_name <- reference_data$tip.label[n]
	tip_name_len <- nchar(tip_name)
	if (tip_name_len > 48) {
		print(tip_name)
		new_tip_name <- str_sub(tip_name, 1, 48)
		reference_data$tip.label[n] <- new_tip_name
	}
}

write.tree(reference_data, file = "Gene_tree_trim_long_tips.txt")