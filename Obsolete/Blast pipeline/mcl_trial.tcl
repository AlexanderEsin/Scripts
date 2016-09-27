direct=~/desktop/Fungi-Archaea-Bact/Rbbh_all

cd $direct

mcl Master_rbh_weight --abc -I 1.2 -o Master_groups_weight -te 7 -scheme 5 --abc-neg-log10 -abc-tf 'ceil(200)' 

######

#-te corresponds to no. threads used
#-scheme (the higher the number the more refined the analysis)
#-abc-neg-log10 allows to use log10 values as edge weights
#-abc-tf (ceil200) caps the highest weight at 200 when 1E-199 = 199