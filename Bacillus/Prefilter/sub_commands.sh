qsub -e ~/err -o ~/out /home/ade110/Scripts/Bacillus/Prefilter/Main/3.Blast_master.sh
qsub -e ~/err -o ~/out /home/ade110/Scripts/Bacillus/Prefilter/Main/4.Blast_out_rearrange_master.sh
qsub -e ~/err -o ~/out /home/ade110/Scripts/Bacillus/Prefilter/Main/5.RBBH_calculate_master.sh
qsub -e ~/err -o ~/out /home/ade110/Scripts/Bacillus/Prefilter/Main/6.RBBH_collate_master.sh
qsub -e ~/err -o ~/out /home/ade110/Scripts/Bacillus/Prefilter/Main/7.Adjust_RBBH_master.sh
qsub -e ~/err -o ~/out /home/ade110/Scripts/Bacillus/Prefilter/Main/8.MCL_master.sh



qsub -e ~/err -o ~/out /home/ade110/Scripts/Bacillus/BLAST/2.Blast_master.sh
qsub -e ~/err -o ~/out /home/ade110/Scripts/Bacillus/RBBH_MCL/3.RBBH_calculate_master.sh
qsub -e ~/err -o ~/out /home/ade110/Scripts/Bacillus/RBBH_MCL/4.RBBH_collate_master.sh
qsub -e ~/err -o ~/out /home/ade110/Scripts/Bacillus/RBBH_MCL/5.0.Adjust_RBBH_master.sh
qsub -e ~/err -o ~/out /home/ade110/Scripts/Bacillus/RBBH_MCL/5.1.MCL_master.sh

qsub -e ~/err -o ~/out /home/ade110/Scripts/Bacillus/Align_tree/9.Group_alignment_master.sh
qsub -e ~/err -o ~/out /home/ade110/Scripts/Bacillus/Align_tree/10.Group_fastTree_master.sh

qsub -e ~/err -o ~/out -v PENALTY='3' /home/ade110/Scripts/Bacillus/Mowgli/2.Run_mowgli_master.sh
qsub -e ~/err -o ~/out -v PENALTY='4' /home/ade110/Scripts/Bacillus/Mowgli/2.Run_mowgli_master.sh
qsub -e ~/err -o ~/out -v PENALTY='4' /home/ade110/Scripts/Bacillus/Mowgli/2.Run_mowgli_master_short.sh
qsub -e ~/err -o ~/out -v PENALTY='3' /home/ade110/Scripts/Bacillus/Mowgli/2.Run_mowgli_master_ax4.sh


qsub -e ~/err -o ~/out /home/ade110/Scripts/Bacillus/Mowgli/0.countCompleteRecon.R



qsub -e ~/err -o ~/out /home/ade110/Scripts/ax4_unpack.sh


###################################

qsub -e ~/err -o ~/out  /mnt/storage/home/aesin/Scripts/kry_misc.sh


qsub -e ~/err -o ~/out -v PENALTY='3' /mnt/storage/home/aesin/Scripts/Bacillus/Mowgli/2.Run_mowgli_master_krypton.sh

qsub -e ~/err -o ~/out /mnt/storage/home/aesin/Scripts/Bacillus/Mowgli/0.countCompleteRecon_krypton