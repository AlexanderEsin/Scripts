qsub -e ~/err -o ~/out /home/ade110/Scripts/Staph/Prefilter/Main/3.Blast_master.sh
qsub -e ~/err -o ~/out /home/ade110/Scripts/Staph/Prefilter/Main/4.Blast_out_rearrange_master.sh
qsub -e ~/err -o ~/out /home/ade110/Scripts/Staph/Prefilter/Main/5.RBBH_calculate_master.sh
qsub -e ~/err -o ~/out /home/ade110/Scripts/Staph/Prefilter/Main/6.RBBH_collate_master.sh
qsub -e ~/err -o ~/out /home/ade110/Scripts/Staph/Prefilter/Main/7.Adjust_RBBH_master.sh
qsub -e ~/err -o ~/out /home/ade110/Scripts/Staph/Prefilter/Main/8.MCL_master.sh


qsub -e ~/err -o ~/out /home/ade110/Scripts/Staph/BLAST/2.Blast_master.sh
qsub -e ~/err -o ~/out /home/ade110/Scripts/Staph/RBBH_MCL/3.RBBH_calculate_master.sh
qsub -e ~/err -o ~/out -v evalue='1e-100' /home/ade110/Scripts/Staph/RBBH_MCL/4.RBBH_collate_master.sh
qsub -e ~/err -o ~/out -v evalue='1e-150' /home/ade110/Scripts/Staph/RBBH_MCL/4.RBBH_collate_master.sh

qsub -e ~/err -o ~/out -v evalue='1e-10' /home/ade110/Scripts/Staph/RBBH_MCL/5.0.Adjust_RBBH_master.sh
qsub -e ~/err -o ~/out -v evalue='1e-50' /home/ade110/Scripts/Staph/RBBH_MCL/5.0.Adjust_RBBH_master.sh
qsub -e ~/err -o ~/out -v evalue='1e-100' /home/ade110/Scripts/Staph/RBBH_MCL/5.0.Adjust_RBBH_master.sh
qsub -e ~/err -o ~/out -v evalue='1e-150' /home/ade110/Scripts/Staph/RBBH_MCL/5.0.Adjust_RBBH_master.sh

qsub -e ~/err -o ~/out -v evalue='-10' /home/ade110/Scripts/Staph/RBBH_MCL/5.1.MCL_master.sh
qsub -e ~/err -o ~/out -v evalue='-50' /home/ade110/Scripts/Staph/RBBH_MCL/5.1.MCL_master.sh
qsub -e ~/err -o ~/out -v evalue='-100' /home/ade110/Scripts/Staph/RBBH_MCL/5.1.MCL_master.sh
qsub -e ~/err -o ~/out -v evalue='-150' /home/ade110/Scripts/Staph/RBBH_MCL/5.1.MCL_master.sh

qsub -e ~/err -o ~/out /home/ade110/Scripts/ax4_unpack.sh