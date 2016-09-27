#!/usr/bin/Rscript
library(ape)
library(phangorn)
library(geiger)

# Counters #
AG_mono_AG_only=0
AG_mono_No_inside=0
AG_mono_All_inside=0
AG_mono_Inside=0
AG_not_mono=0
AG_not_mono_Geo_mono=0
AG_not_mono_Geo_not_mono=0

# Significance values #
G=0
I=0
nT=0

# Prep data frame #
N <- 10000
DF <- data.frame(folder=rep(NA, N), nT=rep(NA, N), G=rep(NA, N), I=rep(NA, N), p=rep(NA, N), stringsAsFactors=FALSE)
DF_sig <- data.frame(folder=rep(NA, N), nT=rep(NA, N), G=rep(NA, N), I=rep(NA, N), p=rep(NA, N), stringsAsFactors=FALSE)

# Bacillac_mono/Anoxy_geo_only with > 8 geo
geo_interesting_list=list()
# List containing Anoxygeo number when there is no inside group
Number_of_anoxy_geo_tips=NULL


tree_build_dir <- "TREE_BUILD"
final_tree_dir <- "FINAL_TREE"
setwd(tree_build_dir)
#print(getwd())
folders <- dir()
i=1

total <- length(folders)

for (folder in folders) {

  cat("\nEstablishing monophylies and tip distributions ...:", i, "/", length(folders), "===", folder)

  setwd(file.path(tree_build_dir, folder))

  # Here we read in the file and root the tree by midpoint for easier visual analysis #
  mp_tree <- Sys.glob(file.path(getwd(), "*_mp_tree*"))
  #print(mp_tree)
  tree <- read.tree(mp_tree)
  done_mptree <- midpoint(tree)
  write.tree(done_mptree,file=mp_tree)
  # We take the original tree and unroot it for the monophyly tests, we will use this from now on #
  unrooted_tree <- unroot(tree)

  # Extract all the tips from the tip files #
  inside_group_tips <- scan(Sys.glob(file.path(getwd(),"*inside_group*")), character(0), quiet=TRUE)
  # For fun extract all the archaea (if any) and then can highlight potential long distance trees #
  archaea <- grep("\\{Arch\\}",unrooted_tree$tip.label,value=T)

  # Get the Geobacilli and Anoxybacilli tips, combine together for the anoxy_geo clade #
  anoxy <- grep("Anoxybacillus",inside_group_tips,value=T)
  geobacs <- grep("\\{Geobac\\}",unrooted_tree$tip.label,value=T)
  anoxy_geo_tips <- c(geobacs,anoxy)

  all_tips <- unrooted_tree$tip.label

  if (is.monophyletic(unrooted_tree,tips=anoxy_geo_tips, reroot=TRUE) == TRUE) {
    if (length(all_tips) == length(anoxy_geo_tips)) {
      file.copy(mp_tree,file.path(final_tree_dir,"Anoxygeo_mono/Anoxy_geo_only"),overwrite=T)
      AG_mono_AG_only=AG_mono_AG_only+1
    } else {
      if (( length(inside_group_tips) - length(anoxy_geo_tips) ) == 0) {
        file.copy(mp_tree,file.path(final_tree_dir,"Anoxygeo_mono/No_inside_group"),overwrite=T)
        AG_mono_No_inside=AG_mono_No_inside+1

        if (length(anoxy_geo_tips) > 8) {
          geo_interesting_list=c(geo_interesting_list, folder)
        }

        # How many anoxy/geobacillus taxa there are in this family
        anoxy_geo_no=length(anoxy_geo_tips)
        Number_of_anoxy_geo_tips=c(Number_of_anoxy_geo_tips, anoxy_geo_no)


      } else if (length(all_tips) == length(inside_group_tips)) {
        file.copy(mp_tree,file.path(final_tree_dir,"Anoxygeo_mono/All_inside_group"),overwrite=T)
        AG_mono_All_inside=AG_mono_All_inside+1

      } else {
        file.copy(mp_tree,file.path(final_tree_dir,"Anoxygeo_mono/Inside_group"),overwrite=T)
        AG_mono_Inside=AG_mono_Inside+1

        G=length(anoxy_geo_tips)
        I=(length(inside_group_tips) - length(anoxy_geo_tips))
        nT=length(all_tips)

        p=((G / I)^2) * (nT / 50)

        DF[i, ] = c(folder, nT, G, I, p)

        if (G > 7 && p > 1 && nT > 50) {
          DF_sig[i, ] = c(folder, nT, G, I, p)
        }

      }
    }
  } else {
    if (is.monophyletic(unrooted_tree, tips = geobacs, reroot=TRUE) == TRUE) {
      file.copy(mp_tree,file.path(final_tree_dir,"Anoxygeo_not_mono/Geo_mono"),overwrite=T)
      AG_not_mono_Geo_mono=AG_not_mono_Geo_mono+1
    } else {
      file.copy(mp_tree,file.path(final_tree_dir,"Anoxygeo_not_mono/Geo_not_mono"),overwrite=T)
      AG_not_mono_Geo_not_mono=AG_not_mono_Geo_not_mono+1
    }
    AG_not_mono=AG_not_mono+1
  }

  i=i+1
}

### Produce dataframe with the "probability" equation ###
DF<-DF[-which(is.na(DF$folder)),]
DF_sig<-DF_sig[-which(is.na(DF_sig$folder)),]
write.table(DF, file="/users/aesin/desktop/Tree_sorting/Test_DF_eqn2.tsv", sep="\t", row.names=F, quote=F)
write.table(DF_sig, file="/users/aesin/desktop/Tree_sorting/Test_DF_eqn2_sig.tsv", sep="\t", row.names=F, quote=F)

### Total trees done ###
cat("\n\nNumber of trees rooted by midpoint:", i-1, "\n")

### Find the numbers of anoxy/geobacillus taxa in each family when there is no inside group ###
table_anoxy_geo_no=as.data.frame(table(Number_of_anoxy_geo_tips))
cat("\nFor those trees with no Inside group the distribution of Anoxy/Geobacillus tips are as follows:\n")
print(table_anoxy_geo_no, row.names=FALSE)


### Counter output ###
cat("\nAG_mono_AG_only:", AG_mono_AG_only, "\nAG_mono_No_inside:", AG_mono_No_inside, "\nAG_mono_All_inside:", AG_mono_All_inside, "\nAG_mono_Inside:", AG_mono_Inside, "\nAG_not_mono:", AG_not_mono, "\n\tAG_not_mono_Geo_mono:", AG_not_mono_Geo_mono, "\n\tAG_not_mono_Geo_not_mono:", AG_not_mono_Geo_not_mono, "\n")
cat("\nNo_inside with > 8 AG:") 
cat(unlist(geo_interesting_list), "\n\n")


# cat("\nNumber of trees with Bacillaceae only:", bacillac_only_count, "\n\tOf which:\n\t\t", bo_anoxy_geo_only_count, "only contained the Anoxy_geobacillus clade\n\t\t", bo_anoxy_geo_mono,"were monophyletic for the Anoxy_geobacillus clade\n\t\t", bo_anoxy_geo_not_mono,"were not monophyletic for the Anoxy_geobacillus clade")

# cat("\nNumber of trees with monophyletic Bacillaceae:", bacillac_mono_count, "\n\tOf which:\n\t\t", bm_all_anoxygeo, "in which n(Bacillaceae) = n(Anoxy_geo)\n\t\t", bm_not_all_anoxygeo, "in which n(Bacillaceae) =/= n(Anoxy_geo)", "\n\t\t\tOf which:\n\t\t\t\t", bm_not_all_anoxygeo_all_bacillales, "in which n(Bacillaceae) = n(Bacillales)\n\t\t\t\t", bm_not_all_anoxygeo_not_all_bacillales, "in which n(Bacillaceae) =/= n(Bacillales)")

# cat("\nNumber of trees with non-monophyletic Bacillaceae:", bacillac_not_mono_count, "\n\tOf which:\n\t\t", bn_anoxy_geo_mono, "were monophyletic for the Anoxy_geobacillus clade", "\n\t\t\tOf which:\n\t\t\t\t", bn_anoxy_geo_mono_bacillal_mono, "were monophyletic for the Bacillales clade\n\t\t\t\t", bn_anoxy_geo_mono_bacillal_not_mono, "were not monophyletic for the Bacillales clade\n\t\t", bn_anoxy_geo_not_mono, "were not monophyletic for the Anoxy_geobacillus clade\n")

# cat(unlist(geo_interesting_list))