library(ape)
library(phangorn)
library(geiger)
library(stringr)

setwd("/users/aesin/desktop/Tree_sorting/Final_trees/Analysis/Bacillac_mono/All_bacillac_anoxygeo")

folders<-dir()
total<-length(folders)

tree_file<-folders[1]
tree<-read.tree(tree_file)
unrooted_tree<-unroot(tree)

all_tips<-unrooted_tree$tip.label

tip_names=list()
for (tip in all_tips) {
	tip_name<-toString(str_sub(str_split_fixed(tip, "\\{", n=2)[1,1],0,-2))
	tip_names<-c(tip_names,tip_name)
}

unlist(tip_names)