cat *result.random_* > all_trees
raxml -f G -m GTRGAMMA -s HBB_no_ovtex.phy -m GTRGAMMA -z all_trees -n logLikes



## Create 100 ML tree forest with constraint
raxml -s HBB_no_ovtex.phy \
-n constrained_rand_step.txt \
-w /users/aesin/Desktop/Deer/Sequences/HBB_vs_HBG/Topology_testing/HBB/Constrained/Random_step/ \
-g Old_new_world_const_rename.txt \
-m GTRGAMMA \
-p $RANDOM \
-N 100 \
-T 18

## Create a good consensus ML tree with constraint
raxml -s HBB_no_ovtex.phy \
-f a \
-n constrained_ML.txt \
-w /users/aesin/Desktop/Deer/Sequences/HBB_vs_HBG/Topology_testing/HBB/Constrained/ML/ \
-g Old_new_world_const_rename.txt \
-m GTRGAMMA \
-p $RANDOM \
-x $RANDOM \
-N 1000 \
-T 18

## Concatenate the true HBBA ML, the constrained HBBA (best) ML, and the 100 ML forest (constrained). Tree 1 = true ML tree; Tree 2 = Best constrained ML tree; Trees 3-102 = constrained forest

cat /users/aesin/Desktop/Deer/Sequences/HBB_vs_HBG/Topology_testing/HBB/ML_tree/*bestTree* /users/aesin/Desktop/Deer/Sequences/HBB_vs_HBG/Topology_testing/HBB/Constrained/ML/*bestTree* /users/aesin/Desktop/Deer/Sequences/HBB_vs_HBG/Topology_testing/HBB/Constrained/Random_step/*result.constrained_rand_step* > /users/aesin/Desktop/Deer/Sequences/HBB_vs_HBG/Topology_testing/HBB/Constrained/all_trees_constrain_test.txt

# Copy the phy alignment file into the testing directory and optimise branch lengths
cd /users/aesin/Desktop/Deer/Sequences/HBB_vs_HBG/Topology_testing/HBB/Constrained
cp /users/aesin/Desktop/Deer/Sequences/HBB_vs_HBG/Topology_testing/HBB/HBB_no_ovtex.phy ./
raxml -f G -s HBB_no_ovtex.phy -m GTRGAMMA -z all_trees_constrain_test.txt -n sitelh


## Run CONSEL
makermt --puzzle RAxML_perSiteLLs
consel RAxML_perSiteLLs
catpv RAxML_perSiteLLs > results.txt