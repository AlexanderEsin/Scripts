#!/usr/bin/Rscript
library(ape)
library(phangorn)
library(geiger)

# Bacillaceae only counts #
bacillac_only_count=0
bo_anoxy_geo_only_count=0
bo_anoxy_geo_mono=0
bo_anoxy_geo_not_mono=0
# Bacillaceae mono counts #
bacillac_mono_count=0
bm_all_anoxygeo=0
bm_not_all_anoxygeo=0
bm_not_all_anoxygeo_all_bacillales=0
bm_not_all_anoxygeo_not_all_bacillales=0
# Bacillaceae not_mono counts #
bacillac_not_mono_count=0
bn_anoxy_geo_mono=0
bn_anoxy_geo_mono_bacillal_mono=0
bn_anoxy_geo_mono_bacillal_not_mono=0
bn_anoxy_geo_not_mono=0
# Bacillac_mono/Anoxy_geo_only with > 8 geo
geo_interesting_list=list()

tree_build_dir<-"TREE_BUILD"
final_tree_dir<-"FINAL_TREE"
setwd(tree_build_dir)
#print(getwd())
folders<-dir()
i=1

total <- length(folders)
# create progress bar
# pb <- tkProgressBar(title = "Calculation progress", min = 0, max = total, width = 300)

for (folder in folders) {
  # setTkProgressBar(pb, i, label=paste( round(i/total*100, 0),"% done"))

  setwd(file.path(tree_build_dir, folder))

  # Here we read in the file and root the tree by midpoint for easier visual analysis #
  mp_tree<-Sys.glob(file.path(getwd(), "*_mp_tree*"))
  #print(mp_tree)
  tree<-read.tree(mp_tree)
  done_mptree<-midpoint(tree)
  write.tree(done_mptree,file=mp_tree)
  # We take the original tree and unroot it for the monophyly tests, we will use this from now on #
  unrooted_tree<-unroot(tree)

  # Extract all the tips from the tip files #
  bacillaceae_tips<-scan(Sys.glob(file.path(getwd(),"*Bacillaceae*")), character(0), quiet=TRUE)
  bacillales_tips<-scan(Sys.glob(file.path(getwd(),"*Bacillales*")), character(0), quiet=TRUE)
  bacilli_tips<-scan(Sys.glob(file.path(getwd(),"*Bacilli*")), character(0), quiet=TRUE)
  firmicutes_tips<-scan(Sys.glob(file.path(getwd(),"*Firmicutes*")), character(0), quiet=TRUE)
  # For fun extract all the archaea (if any) and then can highlight potential long distance trees #
  archaea<-grep("\\{Arch\\}",unrooted_tree$tip.label,value=T)

  # Get the Geobacilli and Anoxybacilli tips, combine together for the anoxy_geo clade #
  anoxy<-grep("Anoxybacillus",bacillaceae_tips,value=T)
  geobacs<-grep("\\{Geobac\\}",unrooted_tree$tip.label,value=T)
  anoxy_geo_tips<-c(geobacs,anoxy)

  all_tips<-unrooted_tree$tip.label

  # If the total number of tips == number of Bacillaceae tips, then all tips are Bacillaceae #
  if (length(all_tips) == length(bacillaceae_tips)) {
    # Within that, if number of all tips == number of Anoxy_geo tips, then all tips are Anoxy_geo #
    if (length(all_tips) == length(anoxy_geo_tips)) {
      file.copy(mp_tree,file.path(final_tree_dir,"Bacillac_only/Geo_Anoxy_only"),overwrite=T)
      bo_anoxy_geo_only_count=bo_anoxy_geo_only_count+1
    # Otherwise check whether the anoxy_geo clade within the Bacillaceae-only tree is monophyletic #
    } else {
      if (is.monophyletic(unrooted_tree,tips=anoxy_geo_tips, reroot=TRUE) == TRUE) {
        file.copy(mp_tree,file.path(final_tree_dir,"Bacillac_only/Geo_Anoxy_mono"),overwrite=T)
        bo_anoxy_geo_mono=bo_anoxy_geo_mono+1
      } else {
        file.copy(mp_tree,file.path(final_tree_dir,"Bacillac_only/Geo_Anoxy_not_mono"),overwrite=T)
        bo_anoxy_geo_not_mono=bo_anoxy_geo_not_mono+1
      }
    }
    #cat("\nDone:", i, "/", length(folders))
    i=i+1
    bacillac_only_count=bacillac_only_count+1
    next
  }

  # If not all the tips are Bacillaceae, check whether the Bacillaceae are monophyletic #
  bacillac_mono<-is.monophyletic(unrooted_tree,tips=bacillaceae_tips, reroot=TRUE)

  # The Bacillaceae are monophyletic #
  if (bacillac_mono == TRUE) {
    if (length(bacillaceae_tips) <= length(anoxy_geo_tips)*1.5) {
      file.copy(mp_tree,file.path(final_tree_dir,"Bacillac_mono/All_bacillac_anoxygeo"),overwrite=T)
      bm_all_anoxygeo=bm_all_anoxygeo+1
      if (length(archaea) > 1) {
        cat("\nArchaea present in Bacillac_mono/All_bacillac_anoxygeo", folder, "no. archaea", length(archaea))
      }
      if (length(anoxy_geo_tips) > 12) {
        geo_interesting_list=c(geo_interesting_list, folder)
      }
    } else {
      if (length(bacillaceae_tips) == length(bacillales_tips)) {
        file.copy(mp_tree,file.path(final_tree_dir,"Bacillac_mono/Not_all_bacillac_anoxygeo/All_bacillac_bacillales"),overwrite=T)
        bm_not_all_anoxygeo_all_bacillales=bm_not_all_anoxygeo_all_bacillales+1
        if (length(archaea) > 1) {
          cat("\nArchaea present in Bacillac_mono/Not_all_bacillac_anoxygeo/All_bacillac_bacillales", folder, "no. archaea", length(archaea))
        }
      } else {
        file.copy(mp_tree,file.path(final_tree_dir,"Bacillac_mono/Not_all_bacillac_anoxygeo/Not_all_bacillac_bacillales"),overwrite=T)
        bm_not_all_anoxygeo_not_all_bacillales=bm_not_all_anoxygeo_not_all_bacillales+1
      }
      bm_not_all_anoxygeo=bm_not_all_anoxygeo+1
    }
    bacillac_mono_count=bacillac_mono_count+1
  # The Bacillaceae are not monophyletic #
  } else if (bacillac_mono == FALSE) {
    if (is.monophyletic(unrooted_tree,tips=anoxy_geo_tips, reroot=TRUE) == TRUE) {
        if (is.monophyletic(unrooted_tree,tips=bacillales_tips, reroot=TRUE) == TRUE) {
          file.copy(mp_tree,file.path(final_tree_dir,"Bacillac_not_mono/Geo_Anoxy_mono/Bacillales_mono"),overwrite=T)
          bn_anoxy_geo_mono_bacillal_mono=bn_anoxy_geo_mono_bacillal_mono+1
        } else {
          file.copy(mp_tree,file.path(final_tree_dir,"Bacillac_not_mono/Geo_Anoxy_mono/Bacillales_not_mono"),overwrite=T)
          bn_anoxy_geo_mono_bacillal_not_mono=bn_anoxy_geo_mono_bacillal_not_mono+1
        }
      bn_anoxy_geo_mono=bn_anoxy_geo_mono+1
      # if (length(archaea) > 0) {
      #   cat("\nArchaea present in Bacillac_not_mono/Geo_Anoxy_mono", folder)
      # }
    } else {
      file.copy(mp_tree,file.path(final_tree_dir,"Bacillac_not_mono/Geo_Anoxy_not_mono"),overwrite=T)
      bn_anoxy_geo_not_mono=bn_anoxy_geo_not_mono+1
      # if (length(archaea) > 0) {
      #   cat("\nArchaea present in Bacillac_not_mono/Geo_Anoxy_not_mono", folder)
      # }
    }
    bacillac_not_mono_count=bacillac_not_mono_count+1
  }

  cat("\nDone:", i, "/", length(folders))
  i=i+1
}

cat("\nNumber of trees rooted by midpoint:", i-1)

cat("\nNumber of trees with Bacillaceae only:", bacillac_only_count, "\n\tOf which:\n\t\t", bo_anoxy_geo_only_count, "only contained the Anoxy_geobacillus clade\n\t\t", bo_anoxy_geo_mono,"were monophyletic for the Anoxy_geobacillus clade\n\t\t", bo_anoxy_geo_not_mono,"were not monophyletic for the Anoxy_geobacillus clade")

cat("\nNumber of trees with monophyletic Bacillaceae:", bacillac_mono_count, "\n\tOf which:\n\t\t", bm_all_anoxygeo, "in which n(Bacillaceae) = n(Anoxy_geo)\n\t\t", bm_not_all_anoxygeo, "in which n(Bacillaceae) =/= n(Anoxy_geo)", "\n\t\t\tOf which:\n\t\t\t\t", bm_not_all_anoxygeo_all_bacillales, "in which n(Bacillaceae) = n(Bacillales)\n\t\t\t\t", bm_not_all_anoxygeo_not_all_bacillales, "in which n(Bacillaceae) =/= n(Bacillales)")

cat("\nNumber of trees with non-monophyletic Bacillaceae:", bacillac_not_mono_count, "\n\tOf which:\n\t\t", bn_anoxy_geo_mono, "were monophyletic for the Anoxy_geobacillus clade", "\n\t\t\tOf which:\n\t\t\t\t", bn_anoxy_geo_mono_bacillal_mono, "were monophyletic for the Bacillales clade\n\t\t\t\t", bn_anoxy_geo_mono_bacillal_not_mono, "were not monophyletic for the Bacillales clade\n\t\t", bn_anoxy_geo_not_mono, "were not monophyletic for the Anoxy_geobacillus clade\n")

cat(unlist(geo_interesting_list))
# close(pb)