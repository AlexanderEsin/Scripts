direct=/scratch/ade110/Geo_v_all

out1=Out
out2=Out_rev

###########################################################################
cd $direct
mkdir -p Out_all

cd $direct/$out1
out_files=?*.tsv
for file in $out_files
do
	mv $file $direct/Out_all
done

cd $direct/$out2
out_files=?*.tsv
for file in $out_files
do
	mv $file $direct/Out_all
done

cd $direct

gzip -r Out_all
cp -r Out_all /csc/rawdata/Warnecke/ADE
gunzip -r Out_all
