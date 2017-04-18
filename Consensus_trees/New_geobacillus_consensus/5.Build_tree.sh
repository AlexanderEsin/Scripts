#!/bin/sh
cd "/Users/aesin/Desktop/Geo_analysis/New_geobacillus_consensus/"
mkdir -p "Tree"
raxml -f a -s Master_concat_alignment.phy -w "/Users/aesin/Desktop/Geo_analysis/New_geobacillus_consensus/Tree" -n concat_geo_tree.txt -m PROTCATAUTO -p $RANDOM -x $RANDOM -N 100 -T 20