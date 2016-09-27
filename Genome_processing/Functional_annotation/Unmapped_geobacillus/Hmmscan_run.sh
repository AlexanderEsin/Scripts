direct=/users/aesin/desktop/Geo_analysis
db_path=$direct/EggNOG/bactNOG_db
input_path=$direct/Geo_ortholog_nucl/Geobac_proteins_no_coords/Fastas
output_path=$direct/Geo_ortholog_nucl/Geobac_proteins_no_coords/Functional_annotation

hmmscan --cpu 6 --tblout $output_path/All_unmapped_bactNOG_hmmscan_out.txt --noali $db_path/bactNOG.txt $input_path/Geobac_not_mapped_fasta.faa
hmmscan --cpu 6 --tblout $output_path/No_part_hyp_bactNOG_hmmscan_out.txt --noali $db_path/bactNOG.txt $input_path/Non_partial_hypothetical_not_mapped.faa