## HBBA AU test - exons only. The exon dataset processed contains only 1 allele of each species (i.e. no O. v. texanus)

# Working directory
master_dir="/Users/aesin/Desktop/Deer/Sequences/HBB_vs_HBG/Topology_testing/HBB/Exon_only"
partition_dir="/Single_partition"

# Align raw fasta sequences into phylip format
mkdir -p $master_dir/Alignment
muscle -in $master_dir/HBB_exons_no_ovtex_raw.fasta -phyiout $master_dir/Alignment/HBB_exons_no_ovtex_aln.phy

## Build an unconstrained ML tree with 1000 BS repeats + optimise the branch lengths / parameters
mkdir -p $master_dir/$partition_dir/Unconstrained_ML

raxml -s $master_dir/Alignment/HBB_exons_no_ovtex_aln.phy \
-f a \
-n unconstrained_ML.txt \
-w $master_dir/$partition_dir/Unconstrained_ML \
-m GTRGAMMA \
-p $RANDOM \
-x $RANDOM \
-N 1000 \
-T 18

# Optimise the branch lengths / parameters of the unconstrained tree
raxml -s $master_dir/Alignment/HBB_exons_no_ovtex_aln.phy \
-f e \
-t $master_dir/$partition_dir/Unconstrained_ML/RAxML_bestTree.unconstrained_ML.txt \
-m GTRGAMMA \
-w $master_dir/$partition_dir/Unconstrained_ML \
-p $RANDOM \
-n unconstrained_optimal.txt \
-T 18

# Build a constrained ML tree with 1000 BS repeats. The tree is constrained to the New World / Old World deer split. Optimise the branch lengths / parameters
mkdir -p $master_dir/$partition_dir/Constrained_ML

raxml -s $master_dir/Alignment/HBB_exons_no_ovtex_aln.phy \
-f a \
-n constrained_ML.txt \
-w $master_dir/$partition_dir/Constrained_ML/ \
-g $master_dir/../../Old_new_world_const_rename.txt \
-m GTRGAMMA \
-p $RANDOM \
-x $RANDOM \
-N 1000 \
-T 18

# Optimise the branch lengths / parameters of the constrained tree
raxml -s $master_dir/Alignment/HBB_exons_no_ovtex_aln.phy \
-f e \
-t $master_dir/$partition_dir/Constrained_ML/RAxML_bestTree.constrained_ML.txt \
-m GTRGAMMA \
-w $master_dir/$partition_dir/Constrained_ML \
-p $RANDOM \
-n constrained_optimal.txt \
-T 18

# Concatenate the true HBBA ML, the constrained HBBA (best) ML, and the 100 ML forest (constrained).
# Tree 1 = true ML tree; Tree 2 = Best constrained ML tree; Trees 3-102 = constrained forest
cat $master_dir/$partition_dir/Unconstrained_ML/RAxML_result.unconstrained_optimal.txt $master_dir/$partition_dir/Constrained_ML/RAxML_result.constrained_optimal.txt > $master_dir/$partition_dir/Trees_concat.txt

# Extract logLikelihood values
mkdir -p $master_dir/$partition_dir/Consel

raxml -s $master_dir/Alignment/HBB_exons_no_ovtex_aln.phy \
-f G \
-w $master_dir/$partition_dir/Consel \
-m GTRGAMMA \
-z $master_dir/$partition_dir/Trees_concat.txt \
-n sitelh \
-T 18

# Run CONSEL
makermt --puzzle $master_dir/$partition_dir/Consel/RAxML_perSiteLLs
consel $master_dir/$partition_dir/Consel/RAxML_perSiteLLs
catpv $master_dir/$partition_dir/Consel/RAxML_perSiteLLs > $master_dir/$partition_dir/Consel/Consel_results.txt


######

partition_dir="/ByBase_partition"

## Build an unconstrained ML tree with 1000 BS repeats + optimise the branch lengths / parameters
mkdir -p $master_dir/$partition_dir/Unconstrained_ML

raxml -s $master_dir/Alignment/HBB_exons_no_ovtex_aln.phy \
-f a \
-n unconstrained_ML.txt \
-w $master_dir/$partition_dir/Unconstrained_ML \
-q $master_dir/byBasePartition.txt \
-m GTRGAMMA \
-p $RANDOM \
-x $RANDOM \
-N 1000 \
-T 18

# Optimise the branch lengths / parameters of the unconstrained tree
raxml -s $master_dir/Alignment/HBB_exons_no_ovtex_aln.phy \
-f e \
-t $master_dir/$partition_dir/Unconstrained_ML/RAxML_bestTree.unconstrained_ML.txt \
-m GTRGAMMA \
-w $master_dir/$partition_dir/Unconstrained_ML \
-q $master_dir/byBasePartition.txt \
-p $RANDOM \
-n unconstrained_optimal.txt \
-T 18

# Build a constrained ML tree with 1000 BS repeats. The tree is constrained to the New World / Old World deer split. Optimise the branch lengths / parameters
mkdir -p $master_dir/$partition_dir/Constrained_ML

raxml -s $master_dir/Alignment/HBB_exons_no_ovtex_aln.phy \
-f a \
-n constrained_ML.txt \
-w $master_dir/$partition_dir/Constrained_ML/ \
-q $master_dir/byBasePartition.txt \
-g $master_dir/../../Old_new_world_const_rename.txt \
-m GTRGAMMA \
-p $RANDOM \
-x $RANDOM \
-N 1000 \
-T 18

# Optimise the branch lengths / parameters of the constrained tree
raxml -s $master_dir/Alignment/HBB_exons_no_ovtex_aln.phy \
-f e \
-t $master_dir/$partition_dir/Constrained_ML/RAxML_bestTree.constrained_ML.txt \
-m GTRGAMMA \
-w $master_dir/$partition_dir/Constrained_ML \
-q $master_dir/byBasePartition.txt \
-p $RANDOM \
-n constrained_optimal.txt \
-T 18

# Concatenate the true HBBA ML, the constrained HBBA (best) ML, and the 100 ML forest (constrained).
# Tree 1 = true ML tree; Tree 2 = Best constrained ML tree; Trees 3-102 = constrained forest
cat $master_dir/$partition_dir/Unconstrained_ML/RAxML_result.unconstrained_optimal.txt $master_dir/$partition_dir/Constrained_ML/RAxML_result.constrained_optimal.txt > $master_dir/$partition_dir/Trees_concat.txt

# Extract logLikelihood values
mkdir -p $master_dir/$partition_dir/Consel

raxml -s $master_dir/Alignment/HBB_exons_no_ovtex_aln.phy \
-f G \
-w $master_dir/$partition_dir/Consel \
-m GTRGAMMA \
-q $master_dir/byBasePartition.txt \
-z $master_dir/$partition_dir/Trees_concat.txt \
-n sitelh \
-T 18

# Run CONSEL
makermt --puzzle $master_dir/$partition_dir/Consel/RAxML_perSiteLLs
consel $master_dir/$partition_dir/Consel/RAxML_perSiteLLs
catpv $master_dir/$partition_dir/Consel/RAxML_perSiteLLs > $master_dir/$partition_dir/Consel/Consel_results.txt


######


partition_dir="/Include_Ovtex"

# Align raw fasta sequences into phylip format
mkdir -p $master_dir/Alignment
muscle -in $master_dir/HBB_exons_raw.fasta -phyiout $master_dir/Alignment/HBB_exons_raw_aln.phy

## Build an unconstrained ML tree with 1000 BS repeats + optimise the branch lengths / parameters
mkdir -p $master_dir/$partition_dir/Unconstrained_ML

raxml -s $master_dir/Alignment/HBB_exons_raw_aln.phy \
-f a \
-n unconstrained_ML.txt \
-w $master_dir/$partition_dir/Unconstrained_ML \
-m GTRGAMMA \
-p $RANDOM \
-x $RANDOM \
-N 1000 \
-T 18

# Optimise the branch lengths / parameters of the unconstrained tree
raxml -s $master_dir/Alignment/HBB_exons_raw_aln.phy \
-f e \
-t $master_dir/$partition_dir/Unconstrained_ML/RAxML_bestTree.unconstrained_ML.txt \
-m GTRGAMMA \
-w $master_dir/$partition_dir/Unconstrained_ML \
-p $RANDOM \
-n unconstrained_optimal.txt \
-T 18


# Build a constrained ML tree with 1000 BS repeats. The tree is constrained to the New World / Old World deer split. Optimise the branch lengths / parameters
mkdir -p $master_dir/$partition_dir/Constrained_ML

raxml -s $master_dir/Alignment/HBB_exons_raw_aln.phy \
-f a \
-n constrained_ML.txt \
-w $master_dir/$partition_dir/Constrained_ML/ \
-g $master_dir/../../Old_new_world_const_ovtex.txt \
-m GTRGAMMA \
-p $RANDOM \
-x $RANDOM \
-N 1000 \
-T 18

# Optimise the branch lengths / parameters of the constrained tree
raxml -s $master_dir/Alignment/HBB_exons_raw_aln.phy \
-f e \
-t $master_dir/$partition_dir/Constrained_ML/RAxML_bestTree.constrained_ML.txt \
-m GTRGAMMA \
-w $master_dir/$partition_dir/Constrained_ML \
-p $RANDOM \
-n constrained_optimal.txt \
-T 18

# Concatenate the true HBBA ML, the constrained HBBA (best) ML, and the 100 ML forest (constrained).
# Tree 1 = true ML tree; Tree 2 = Best constrained ML tree; Trees 3-102 = constrained forest
cat $master_dir/$partition_dir/Unconstrained_ML/RAxML_result.unconstrained_optimal.txt $master_dir/$partition_dir/Constrained_ML/RAxML_result.constrained_optimal.txt > $master_dir/$partition_dir/Trees_concat.txt

# Extract logLikelihood values
mkdir -p $master_dir/$partition_dir/Consel

raxml -s $master_dir/Alignment/HBB_exons_raw_aln.phy \
-f G \
-w $master_dir/$partition_dir/Consel \
-m GTRGAMMA \
-z $master_dir/$partition_dir/Trees_concat.txt \
-n sitelh \
-T 18

# Run CONSEL
makermt --puzzle $master_dir/$partition_dir/Consel/RAxML_perSiteLLs
consel $master_dir/$partition_dir/Consel/RAxML_perSiteLLs
catpv $master_dir/$partition_dir/Consel/RAxML_perSiteLLs > $master_dir/$partition_dir/Consel/Consel_results.txt


######

partition_dir="/Pig_outgroup"

# Align raw fasta sequences into phylip format
mkdir -p $master_dir/Alignment
muscle -in $master_dir/HBB_exons_outgroup_raw.fasta -phyiout $master_dir/Alignment/HBB_exons_outgroup_raw.phy

## Build an unconstrained ML tree with 1000 BS repeats + optimise the branch lengths / parameters
mkdir -p $master_dir/$partition_dir/Unconstrained_ML

raxml -s $master_dir/Alignment/HBB_exons_outgroup_raw.phy \
-f a \
-n unconstrained_ML.txt \
-w $master_dir/$partition_dir/Unconstrained_ML \
-o Sscro_HBB \
-m GTRGAMMA \
-p $RANDOM \
-x $RANDOM \
-N 1000 \
-T 18

# Optimise the branch lengths / parameters of the unconstrained tree
raxml -s $master_dir/Alignment/HBB_exons_outgroup_raw.phy \
-f e \
-t $master_dir/$partition_dir/Unconstrained_ML/RAxML_bestTree.unconstrained_ML.txt \
-m GTRGAMMA \
-w $master_dir/$partition_dir/Unconstrained_ML \
-p $RANDOM \
-n unconstrained_optimal.txt \
-T 18


# Build a constrained ML tree with 1000 BS repeats. The tree is constrained to the New World / Old World deer split. Optimise the branch lengths / parameters
mkdir -p $master_dir/$partition_dir/Constrained_ML

raxml -s $master_dir/Alignment/HBB_exons_outgroup_raw.phy \
-f a \
-n constrained_ML.txt \
-w $master_dir/$partition_dir/Constrained_ML/ \
-o Sscro_HBB \
-g $master_dir/../../Old_new_world_const_sscrofa.txt \
-m GTRGAMMA \
-p $RANDOM \
-x $RANDOM \
-N 1000 \
-T 18

# # Optimise the branch lengths / parameters of the constrained tree
raxml -s $master_dir/Alignment/HBB_exons_outgroup_raw.phy \
-f e \
-t $master_dir/$partition_dir/Constrained_ML/RAxML_bestTree.constrained_ML.txt \
-m GTRGAMMA \
-w $master_dir/$partition_dir/Constrained_ML \
-p $RANDOM \
-n constrained_optimal.txt \
-T 18

# Concatenate the true HBBA ML, the constrained HBBA (best) ML, and the 100 ML forest (constrained).
# Tree 1 = true ML tree; Tree 2 = Best constrained ML tree; Trees 3-102 = constrained forest
cat $master_dir/$partition_dir/Unconstrained_ML/RAxML_result.unconstrained_optimal.txt $master_dir/$partition_dir/Constrained_ML/RAxML_result.constrained_optimal.txt > $master_dir/$partition_dir/Trees_concat.txt

# Extract logLikelihood values
mkdir -p $master_dir/$partition_dir/Consel

raxml -s $master_dir/Alignment/HBB_exons_outgroup_raw.phy \
-f G \
-w $master_dir/$partition_dir/Consel \
-m GTRGAMMA \
-z $master_dir/$partition_dir/Trees_concat.txt \
-n sitelh \
-T 18

# Run CONSEL
makermt --puzzle $master_dir/$partition_dir/Consel/RAxML_perSiteLLs
consel $master_dir/$partition_dir/Consel/RAxML_perSiteLLs
catpv $master_dir/$partition_dir/Consel/RAxML_perSiteLLs > $master_dir/$partition_dir/Consel/Consel_results.txt


##########

partition_dir="/New_const"

# Align raw fasta sequences into phylip format
mkdir -p $master_dir/Alignment
muscle -in $master_dir/HBB_exons_raw.fasta -phyiout $master_dir/Alignment/HBB_exons_raw_aln.phy

## Build an unconstrained ML tree with 1000 BS repeats + optimise the branch lengths / parameters
mkdir -p $master_dir/$partition_dir/Unconstrained_ML

raxml -s $master_dir/Alignment/HBB_exons_raw_aln.phy \
-f a \
-n unconstrained_ML.txt \
-w $master_dir/$partition_dir/Unconstrained_ML \
-m GTRGAMMA \
-p $RANDOM \
-x $RANDOM \
-N 1000 \
-T 18

# Optimise the branch lengths / parameters of the unconstrained tree
raxml -s $master_dir/Alignment/HBB_exons_raw_aln.phy \
-f e \
-t $master_dir/$partition_dir/Unconstrained_ML/RAxML_bestTree.unconstrained_ML.txt \
-m GTRGAMMA \
-w $master_dir/$partition_dir/Unconstrained_ML \
-p $RANDOM \
-n unconstrained_optimal.txt \
-T 18

# Build a constrained ML tree with 1000 BS repeats. The tree is constrained to the New World / Old World deer split. Optimise the branch lengths / parameters
mkdir -p $master_dir/$partition_dir/Constrained_ML

raxml -s $master_dir/Alignment/HBB_exons_raw_aln.phy \
-f a \
-n constrained_ML.txt \
-w $master_dir/$partition_dir/Constrained_ML/ \
-g $master_dir/../../NEW_const.txt \
-m GTRGAMMA \
-p $RANDOM \
-x $RANDOM \
-N 1000 \
-T 18

# # Optimise the branch lengths / parameters of the constrained tree
raxml -s $master_dir/Alignment/HBB_exons_raw_aln.phy \
-f e \
-t $master_dir/$partition_dir/Constrained_ML/RAxML_bestTree.constrained_ML.txt \
-m GTRGAMMA \
-w $master_dir/$partition_dir/Constrained_ML \
-p $RANDOM \
-n constrained_optimal.txt \
-T 18

# Concatenate the true HBBA ML, the constrained HBBA (best) ML, and the 100 ML forest (constrained).
# Tree 1 = true ML tree; Tree 2 = Best constrained ML tree; Trees 3-102 = constrained forest
cat $master_dir/$partition_dir/Unconstrained_ML/RAxML_result.unconstrained_optimal.txt $master_dir/$partition_dir/Constrained_ML/RAxML_result.constrained_optimal.txt > $master_dir/$partition_dir/Trees_concat.txt

# Extract logLikelihood values
mkdir -p $master_dir/$partition_dir/Consel

raxml -s $master_dir/Alignment/HBB_exons_raw_aln.phy \
-f G \
-w $master_dir/$partition_dir/Consel \
-m GTRGAMMA \
-z $master_dir/$partition_dir/Trees_concat.txt \
-n sitelh \
-T 18

# Run CONSEL
makermt --puzzle $master_dir/$partition_dir/Consel/RAxML_perSiteLLs
consel $master_dir/$partition_dir/Consel/RAxML_perSiteLLs
catpv $master_dir/$partition_dir/Consel/RAxML_perSiteLLs > $master_dir/$partition_dir/Consel/Consel_results.txt











# # Create 100 ML tree forest using random stepwise addition starting trees. These trees are constrained to the New World / Old World deer split
# mkdir -p $master_dir/Constrained_forest_random

# raxml -s $master_dir/Alignment/HBB_exons_no_ovtex_aln.phy \
# -d \
# -n constrained_rand.txt \
# -w $master_dir/Constrained_forest_random \
# -g $master_dir/../../Old_new_world_const_rename.txt \
# -m GTRGAMMA \
# -p $RANDOM \
# -N 100 \
# -T 18