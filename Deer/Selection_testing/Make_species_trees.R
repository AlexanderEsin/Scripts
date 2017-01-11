library(ape)
library(phangorn)
library(phylobase)
library(geiger)

###########################################################################

setwd("/users/aesin/Desktop/Deer/Selection_analysis/Species_tree/")

timetree <- read.tree("TimetreeOfLife.nwk")
## All tips ##
all_keep_tips <- c("Ailuropoda_melanoleuca", "Bos_taurus", "Camelus_ferus", "Felis_catus", "Myotis_lucifugus", "Homo_sapiens", "Sus_scrofa", "Ovis_aries", "Capra_hircus", "Mus_musculus", "Oryctolagus_cuniculus", "Callithrix_jacchus", "Macaca_mulatta", "Papio_anubis", "Rattus_norvegicus", "Pteropus_vampyrus", "Odobenus_rosmarus", "Odocoileus_virginianus", "Alces_alces", "Cervus_elaphus", "Elaphurus_davidianus", "Cervus_nippon", "Rangifer_tarandus", "Dama_dama", "Moschus_berezovskii", "Equus_caballus", "Ceratotherium_simum", "Eptesicus_fuscus", "Panthera_tigris", "Capreolus_capreolus", "Pudu_puda", "Giraffa_camelopardalis", "Hydropotoes_inermis", "Cervus_albirostris", "Rucervus_duvaucelii", "Muntiacus_reevesi", "Vicugna_pacos")
#"Dasypus_novemcinctus"

for_dating_tree <- drop.tip(timetree, setdiff(timetree$tip.label, all_keep_tips))

##################################

## No deer ##
keep_tips <- c("Ailuropoda_melanoleuca", "Bos_taurus", "Camelus_ferus","Vicugna_pacos", "Felis_catus", "Myotis_lucifugus", "Homo_sapiens", "Sus_scrofa", "Ovis_aries", "Capra_hircus", "Mus_musculus", "Oryctolagus_cuniculus", "Callithrix_jacchus", "Macaca_mulatta", "Papio_anubis", "Rattus_norvegicus", "Pteropus_vampyrus", "Odobenus_rosmarus", "Moschus_berezovskii", "Equus_caballus", "Ceratotherium_simum", "Eptesicus_fuscus", "Panthera_tigris", "Giraffa_camelopardalis", "Vicugna_pacos")
#"Dasypus_novemcinctus"

time_prune_tree <- drop.tip(timetree, setdiff(timetree$tip.label, keep_tips))

time_prune_tree$edge.length[which(time_prune_tree$edge.length == 0)] <- mean(time_prune_tree$edge.length)

##################################

## Change position of carnivores to match Meredith (2011) ##
# Isolate carnivore subtree #
carnivore_node <- getMRCA(time_prune_tree, c("Felis_catus", "Ailuropoda_melanoleuca"))
carnivore_tips <- tips(time_prune_tree, carnivore_node)
carnivore_subtree <- drop.tip(time_prune_tree, setdiff(time_prune_tree$tip.label, carnivore_tips))
carnivore_topo <- chronopl(carnivore_subtree, 0, age.min = 1, age.max = NULL,node = "root", S = 1, tol = 1e-8,CV = FALSE, eval.max = 500, iter.max = 500)

time_prune_tree <- drop.tip(time_prune_tree, carnivore_tips)

		
# Rebind subtree #
message("Bind Carnivora")
plot.phylo(time_prune_tree)
bind_where <- getMRCA(time_prune_tree, c("Myotis_lucifugus", "Ceratotherium_simum"))
rebind_carnivore <- bind.tree(time_prune_tree, carnivore_topo, where = bind_where, position = 0.1)	
rebind_carnivore <- multi2di(rebind_carnivore, random = FALSE)
time_prune_tree <- rebind_carnivore
#time_prune_tree$edge.length <- rep(1, length(time_prune_tree$edge.length))
# rebind_carnivore$edge.length[which(rebind_carnivore$edge.length == 0)] <- 1
#rebind_carnivore <- chronopl(rebind_carnivore, 0, age.min = 1, age.max = NULL,node = "root", S = 1, tol = 1e-8,CV = FALSE, eval.max = 500, iter.max = 500)

##################################

## Get the tenk tree ##
tenk_tree <- read.nexus(file = "Artiodactyl_with_deer_10ktrees.nex")

## Extract deer subtree ##
deer_node <- getMRCA(tenk_tree, c("Dama_dama", "Alces_alces"))
deer_tips <- tips(tenk_tree, deer_node)
deer_prune <- drop.tip(tenk_tree, setdiff(tenk_tree$tip.label, deer_tips))

## Add the Bactrian tip ##
bactrian_tip <- list(edge=matrix(c(2,1),1,2), tip.label="Cervus_elaphus_bactrianus", edge.length=1.0, Nnode=1)
class(bactrian_tip)<-"phylo"
message("Bind the Bactrian deer")
plot.phylo(deer_prune)
deer_prune <- bind.tree(deer_prune, bactrian_tip, interactive = TRUE)
deer_prune <- chronopl(deer_prune, 0, age.min = 1, age.max = NULL,node = "root", S = 1, tol = 1e-8,CV = FALSE, eval.max = 500, iter.max = 500)

## Add the Cervus_canadensis tip ##
wapiti_tip <- list(edge=matrix(c(2,1),1,2), tip.label="Cervus_canadensis", edge.length=1.0, Nnode=1)
class(wapiti_tip)<-"phylo"
message("Bind the Waptiti (sister to nippon")
plot.phylo(deer_prune)
deer_prune <- bind.tree(deer_prune, wapiti_tip, interactive = TRUE)
deer_prune <- chronopl(deer_prune, 0, age.min = 1, age.max = NULL,node = "root", S = 1, tol = 1e-8,CV = FALSE, eval.max = 500, iter.max = 500)

## Rename the Swamp_deer to Rucervus ##
swamp_deer_index <- grep("Cervus_duvaucelii", deer_prune)
deer_prune$tip.label[swamp_deer_index] <- "Rucervus_duvaucelii"

## Ultrametric deer ## 
deer_topo <- chronopl(deer_prune, 0, age.min = 1, age.max = NULL,node = "root", S = 1, tol = 1e-8,CV = FALSE, eval.max = 500, iter.max = 500)
# Write out just the deer tree #
just_deer_write <- deer_topo; just_deer_write$edge.length <- NULL; just_deer_write$node.label <- NULL
write.tree(just_deer_write, file = "Just_deer/Just_deer_topo_tree.txt")

##################################

## Bind deer tree to main tree ##
message("Bind Deer to Artiodactyla")
plot.phylo(time_prune_tree)
bind_where <- getMRCA(time_prune_tree, c("Moschus_berezovskii", "Bos_taurus"))
bound_deer <- bind.tree(time_prune_tree, deer_topo, where = bind_where, position = 0.1)
bound_deer <- multi2di(bound_deer, random = FALSE)
bound_deer$edge.length[which(bound_deer$edge.length == 0)] <- mean(bound_deer$edge.length)

## Set edge legths to NULL to make topo only ##
all_topo <- bound_deer
all_topo$edge.length <- NULL
all_topo$node.label <- NULL


## Write tree ##
write.tree(all_topo, file = "Timetree_prune_meredith.txt")

###########################################################################
## Also write out just the ruminant topology ##

## Extract ruminant subtree ##
rum_node <- getMRCA(all_topo, c("Giraffa_camelopardalis", "Alces_alces"))
rum_tips <- tips(all_topo, rum_node)
rum_prune <- drop.tip(all_topo, setdiff(all_topo$tip.label, rum_tips))

write.tree(rum_prune, file = "Rum_only_topo.txt")

###########################################################################
## Also write out just the deer topology ##

## Extract ruminant subtree ##
deer_node <- getMRCA(all_topo, c("Dama_dama", "Alces_alces"))
deer_tips <- tips(all_topo, deer_node)
deer_prune <- drop.tip(all_topo, setdiff(all_topo$tip.label, deer_tips))

write.tree(deer_prune, file = "Deer_only_topo.txt")

###########################################################################
## Make two reindeer ##

# two_rein_tree_topo <- all_topo
# two_rein_tree_topo$tip.label[which(two_rein_tree_topo$tip.label == "Rangifer_tarandus")] <- "Rangifer_tarandus_fennicus"

# new_rein_tip <- list(edge=matrix(c(2,1),1,2), tip.label="Rangifer_tarandus_tarandus", edge.length=1.0, Nnode=1)
# class(new_rein_tip)<-"phylo"

# # attach to any node (say, node 16)
# plot.phylo(two_rein_tree_topo)
# two_rein_tree_topo <- bind.tree(two_rein_tree_topo, new_rein_tip, interactive = TRUE)
# two_rein_tree_topo <- chronopl(two_rein_tree_topo, 0, age.min = 1, age.max = NULL,node = "root", S = 1, tol = 1e-8,CV = FALSE, eval.max = 500, iter.max = 500)
# plot.phylo(two_rein_tree_topo)

# write.tree(two_rein_tree_topo, "Meredith_topo_two_rein.txt")

 













