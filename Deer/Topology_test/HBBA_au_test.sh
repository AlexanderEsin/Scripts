# Directory : /users/aesin/Desktop/Deer/Sequences/HBB_vs_HBG/Topology_testing/HBB/ML_tree/


# ## Create 100 ML tree forest with constraint
# raxml -s HBB_no_ovtex.phy \
# -n constrained_rand_step.txt \
# -w /users/aesin/Desktop/Deer/Sequences/HBB_vs_HBG/Topology_testing/HBB/Constrained/Random_step/ \
# -g Old_new_world_const_rename.txt \
# -m GTRGAMMA \
# -p $RANDOM \
# -N 100 \
# -T 18

cd /users/aesin/Desktop/Deer/Sequences/HBB_vs_HBG/Topology_testing/HBB/Full_gene/
mkdir -p Unconstrained_ML
mkdir -p Constrained_ML
mkdir -p Consel

## Create a good consensus ML tree with constraint
raxml -s HBB_no_ovtex.phy \
-f a \
-n unconstrained_ML.txt \
-w /users/aesin/Desktop/Deer/Sequences/HBB_vs_HBG/Topology_testing/HBB/Full_gene/Unconstrained_ML/ \
-m GTRGAMMA \
-p $RANDOM \
-x $RANDOM \
-N 1000 \
-T 18

## Create a good consensus ML tree with constraint
raxml -s HBB_no_ovtex.phy \
-f a \
-n constrained_ML.txt \
-w /users/aesin/Desktop/Deer/Sequences/HBB_vs_HBG/Topology_testing/HBB/Full_gene/Constrained_ML/ \
-g Old_new_world_const_rename.txt \
-m GTRGAMMA \
-p $RANDOM \
-x $RANDOM \
-N 1000 \
-T 18

## Concatenate the true HBBA ML, the constrained HBBA (best) ML, and the 100 ML forest (constrained). Tree 1 = true ML tree; Tree 2 = Best constrained ML tree; Trees 3-102 = constrained forest

cat /users/aesin/Desktop/Deer/Sequences/HBB_vs_HBG/Topology_testing/HBB/Full_gene/Unconstrained_ML/*bestTree* /users/aesin/Desktop/Deer/Sequences/HBB_vs_HBG/Topology_testing/HBB/Full_gene/Constrained_ML/*bestTree* > /users/aesin/Desktop/Deer/Sequences/HBB_vs_HBG/Topology_testing/HBB/Full_gene/Consel/combined_trees.txt

# Copy the phy alignment file into the testing directory and optimise branch lengths
cd Consel
raxml -f G -s /users/aesin/Desktop/Deer/Sequences/HBB_vs_HBG/Topology_testing/HBB/Full_gene/HBB_no_ovtex.phy -m GTRGAMMA -z combined_trees.txt -n sitelh


## Run CONSEL
makermt --puzzle RAxML_perSiteLLs
consel RAxML_perSiteLLs
catpv RAxML_perSiteLLs > No_ovtex_HBB_AU_results.txt