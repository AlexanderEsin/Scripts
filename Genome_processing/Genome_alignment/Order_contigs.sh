## Reorder the contigs in the not-fully-assembled genomes to identify those contigs (at the end) that do not align. For reference, the closest fully assembled genome is used. Closeness is based on the Geobacillus consensus species tree ##

cd /users/aesin/desktop/Geo_analysis/Geo_omes/Genome_processing

## G11MC16 ##
# Vs Thermodenitrificans #
mkdir -p ./Contig_reorder_genbank/G11MC16_contig_move_thermodenitrificans
java -Xmx4096m -cp "/Applications/Mauve.app/Contents/Java/*" org.gel.mauve.contigs.ContigOrderer -output ./Contig_reorder_genbank/G11MC16_contig_move_thermodenitrificans -ref ./Geobac_genomes_assembled_genbank/GCF_000015745.1_ASM1574v1_genomic.gbk -draft ./Geobac_genomes_contigs_genbank/GCF_000173035.1_ASM17303v1_genomic.gbk

## Caldoxylosilyticus ##
# Vs WCH70 #
java -Xmx4096m -cp "/Applications/Mauve.app/Contents/Java/*" org.gel.mauve.contigs.ContigOrderer -output ./Contig_reorder_genbank/Caldoxylosilyticus_contig_move_WCH70 -ref ./Geobac_genomes_assembled_genbank/GCF_000023385.1_ASM2338v1_genomic.gbk -draft ./Geobac_genomes_contigs_genbank/GCF_000313345.1_CIC9_01_genomic.gbk

## WSUCF1 ##
# Vs GHH01 #
java -Xmx4096m -cp "/Applications/Mauve.app/Contents/Java/*" org.gel.mauve.contigs.ContigOrderer -output ./Contig_reorder/WSUCF1_contig_move_GHH01 -ref ./Geobac_genomes_assembled_genbank/GCF_000336445.1_ASM33644v1_genomic.gbk -draft ./Geobac_genomes_contigs_genbank/GCF_000422025.1_DUJ_genomic.gbk

## MAS1 ##
# Vs thermoleovorans #
java -Xmx4096m -cp "/Applications/Mauve.app/Contents/Java/*" org.gel.mauve.contigs.ContigOrderer -output ./Contig_reorder/MAS1_contig_move_thermoleovorans -ref ./Geobac_genomes_assembled_genbank/GCF_000236605.1_ASM23660v1_genomic.gbk -draft ./Geobac_genomes_contigs_genbank/GCF_000498995.1_GTMas1.0_genomic.gbk

## Thermocatenulatus ##
# Vs GHH01 #
java -Xmx4096m -cp "/Applications/Mauve.app/Contents/Java/*" org.gel.mauve.contigs.ContigOrderer -output ./Contig_reorder/Thermocatenulatus_contig_move_GHH01 -ref ./Geobac_genomes_assembled_genbank/GCF_000336445.1_ASM33644v1_genomic.gbk -draft ./Geobac_genomes_contigs_genbank/GCF_000612265.1_GeoThe1.0_genomic.gbk

## Stearothermophilus ##
# Vs vulcani #
java -Xmx4096m -cp "/Applications/Mauve.app/Contents/Java/*" org.gel.mauve.contigs.ContigOrderer -output ./Contig_reorder/Stearothermophilus_contig_move_vulcani -ref ./Geobac_genomes_assembled_genbank/GCF_000733845.1_ASM73384v1_genomic.gbk -draft ./Geobac_genomes_contigs_genbank/GCF_000705495.1_GbsDonk1.0_genomic.gbk

## G1W1 ##
# Vs vulcani #
java -Xmx4096m -cp "/Applications/Mauve.app/Contents/Java/*" org.gel.mauve.contigs.ContigOrderer -output ./Contig_reorder/G1W1_contig_move_vulcani -ref ./Geobac_genomes_assembled_genbank/GCF_000733845.1_ASM73384v1_genomic.gbk -draft ./Geobac_genomes_contigs_genbank/GCF_000750005.1_ASM75000v1_genomic.gbk