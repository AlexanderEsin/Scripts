library(ape)
library(phylobase)
library(phytools)
library(stringr)
library(geiger)
library(RSvgDevice)
library(reshape2)
library(ggplot2)
library(gtools)

### 

direct <- "/users/aesin/desktop/Mowgli/"
test_SPR_dir <- paste0(direct, "Test_SPR/")

inside_group_tips	<- as.vector(read.delim(paste0(direct, "Inside_group_tips.txt"), header = FALSE)$V1)
species_tree <- read.tree(file = paste0(direct, "Species_tree/Ultrametric_species_tree_tipfix.txt"))

input_dir	<- paste0(test_SPR_dir, "Input/")
output_dir	<- paste0(test_SPR_dir, "Output/")

#species_tree_dir <- paste0(direct, "Species_tree")
if (file.exists(output_dir) == FALSE) {
	dir.create(output_dir)
}



## Tip names were changed for the species tree to be used by Mowgli. But, the original name is always a subset of the name used in the species tree, so we can find the list of modified labels with grep ##
renamed_IG_tips		<- vector()
for (tip in inside_group_tips) {

	renamed_tip		<- grep(tip, species_tree$tip.label, value = TRUE)
	renamed_IG_tips	<- c(renamed_IG_tips, renamed_tip)
}

setwd(input_dir)
input_trees <- mixedsort(dir())

for (tree_file in input_trees) {
	
	tree_number <- gsub("[^0-9]+", "", tree_file)
	dir.create(paste0(output_dir, tree_number), showWarnings = FALSE)

	print(tree_number)
	full_tree <- read.tree(tree_file)
	print(paste0("Original number of tips: ", length(full_tree$tip.label)))


	IG_present <- vector()
	for (tip in inside_group_tips) {
		hit <- (grep(tip, full_tree$tip.label, value = T))
		if (length(hit) != 0) {
			IG_present <- c(IG_present, hit)
		}
	}

	not_IG_tips <- setdiff(full_tree$tip.label, IG_present)

	for (i in 1:10) {
		out_name <- paste0("spr", i)

		if ((length(IG_present)) >= length(not_IG_tips)) {
			write.tree(full_tree, file = paste0(output_dir, tree_number, "/", out_name, ".txt"))
			if (i == 1) {
				print("All tips kept")
			}
			next
		} else if ((length(IG_present)) < 100) {
			sample_num <- 100
		} else {
			sample_num <- length(IG_present)
		}
	
		random_select <- sample(not_IG_tips, sample_num)
		keep_tips <- c(IG_present, random_select)
		prune <- drop.tip(full_tree, setdiff(full_tree$tip.label, keep_tips))
		write.tree(prune, file = paste0(output_dir, tree_number, "/", out_name, ".txt"))

		if (i == 1) {
			print(paste0("Tips to be reconciled: ", length(keep_tips)))
		}

	}
	print("\n")

}


## 