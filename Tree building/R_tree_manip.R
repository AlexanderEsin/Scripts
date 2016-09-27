rm(list=ls())

library(geiger)
library(ape)
#Pull out any tip_name which has the class {Geobac}#
geobacs<-grep("\\{Geobac\\}",tree$tip.label,value=T)
subtilis<-grep("Bacillus_subtilis",tree$tip.label,value=T)
geobacs_subt<-c(geobacs,subtilis)

arch<-grep("\\{Arch\\}",tree$tip.label,value=T)


#Get the node value for the mrca of all the geobac#
node1=getMRCA(tree,geobacs)
node2=getMRCA(tree,geobacs_subt)

node2=getMRCA(tree,arch)

#Get all the tips subtended by the node#
clade=tips(tree,node2)

#Extract the desired clade based on the node#
bacilli=extract.clade(tree,node2)

#Plot the clade#
pdf("test_tree.pdf")
plot.phylo(bacilli,cex=0.2,no.margin=T)
dev.off()

#We could then compare the geobacs_subt clade with the same clade from the reference tree#
