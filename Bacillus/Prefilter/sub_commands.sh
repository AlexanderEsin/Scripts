qsub -e ~/err -o ~/out /home/ade110/Scripts/Bacillus/Prefilter/Main/3.Blast_master.sh
qsub -e ~/err -o ~/out /home/ade110/Scripts/Bacillus/Prefilter/Main/4.Blast_out_rearrange_master.sh
qsub -e ~/err -o ~/out /home/ade110/Scripts/Bacillus/Prefilter/Main/5.RBBH_calculate_master.sh
qsub -e ~/err -o ~/out /home/ade110/Scripts/Bacillus/Prefilter/Main/6.RBBH_collate_master.sh
qsub -e ~/err -o ~/out /home/ade110/Scripts/Bacillus/Prefilter/Main/7.Adjust_RBBH_master.sh
qsub -e ~/err -o ~/out /home/ade110/Scripts/Bacillus/Prefilter/Main/8.MCL_master.sh



qsub -e ~/err -o ~/out $HOME/Scripts/ax4_unpack.sh