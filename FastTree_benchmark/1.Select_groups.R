#!/usr/local/bin/Rscript

## In this Script we sample 200 lHGT gene families, 100 sHGT families
## and 300 vertical gene families. For each of these gene families
## we will reconstruct the FastTree gene tree and reconcile that gene
## tree with Mowgli. The expectation is that the fastTree based
## reconciliation will show almost identical results as that based
## on the RAxML dataset. This would mean we can then use FastTree instead
## of RAxML for the larger AnoGeo dataset.

## We will take lHGT and sHGT at penalty = 5 (i.e. consistent over 3, 4 ,5)
direct <-		"/Users/aesin/Desktop/Mowgli/"
lHGT_file <-	file.path(direct, "Refined_events", "Events", "T5_full_lHGT_events.tsv")
sHGT_file <-	file.path(direct, "Refined_events", "Events", "T5_full_sHGT_events.tsv")

## Vertical file consistent over 4, 5, 6
vert_file <-	file.path(direct, "Constant_events", "Ver_const_t4_t5_t6.tsv")

## Read in each file
vert_events <-	read.table(vert_file, header = FALSE, sep = "")
lHGT_events <-	read.table(lHGT_file, header = TRUE, sep = "\t")
sHGT_events <-	read.table(sHGT_file, header = TRUE, sep = "\t")

## Sample the requisite number of groups from each type
vert_sample <-	sample(vert_events$V1, 300)
lHGT_sample <-	sample(unique(lHGT_events$Group), 200)
sHGT_sample <-	sample(unique(sHGT_events$Group), 100)

## Combine HGT samples and look for group overlap
HGT_sample	<-	c(lHGT_sample, sHGT_sample)
sHGT_lHGT	<-	HGT_sample[duplicated(HGT_sample)]
# Groups 1276, 1946, 1672 have both short and long HGTs predicted


## A total of 597 unique groups to be processed.
total_redo	<- unique(c(vert_sample, HGT_sample))

write.table(x = total_redo, file = "/Users/aesin/Desktop/FastTree/Input_groups/Group_list.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)


