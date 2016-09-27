#!/bin/sh
##PBS -l walltime=00:30:00
##PBS -l mem=4gb
##PBS -l ncpus=8

##module load blast+/2.2.28

SCRATCH=~/desktop/

cd $SCRATCH/Archaea_v_Fungi/archaea_faa
mkdir -p DB
mkdir -p out
mkdir -p timing

cd $SCRATCH/Archaea_v_Fungi/fungi_db

dblist=?*.faa

for db in $dblist
do
	makeblastdb -in $db -dbtype prot -out $SCRATCH/Archaea_v_Fungi/archaea_faa/DB/$db.db
done


echo "cd $SCRATCH/Archaea_v_Fungi/fungi_db
set dblist [lrange [glob *.faa] [lindex \$argv 0] [lindex \$argv 1]]
cd $SCRATCH/Archaea_v_Fungi/archaea_faa
set glist [glob *.faa]

proc openfile {fl} {
	global data
	set in [open \$fl r]
	set data [read \$in]
	close \$in
	return
}

foreach db \$dblist {
	set err "x\\n"
	foreach g \$glist {
		if {\$g == \$db} {
			puts "skip"
		} else {
			set systime [clock seconds]
			append err "\$g Start: [clock format \$systime -format %H:%M:%S]\\n"
			regsub _protein.faa \$g {} x
			regsub _protein.faa \$db {} y
			catch {exec blastp -query \$g -db DB/\$db.db -out out/\$x\&\$y.tsv -evalue 1e-10 -outfmt 7 -max_target_seqs 1 -max_hsps 1 -num_threads 2}
			set systime [clock seconds]
			append err "\$g Finish: [clock format \$systime -format %H:%M:%S]\\n\\n"
		}
	}
	set out [open $SCRATCH/Archaea_v_Fungi/archaea_faa/timing/\$db.txt a]
	puts \$out \$err
	close \$out
}" > $SCRATCH/arch_v_fung_script.tcl

chmod +x $SCRATCH/arch_v_fung_script.tcl

tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 0 0 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 1 1 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 2 2 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 3 3 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 4 4 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 5 5 &

wait

##rm -r $SCRATCH/Archaea_v_Fungi/archaea_faa/DB
mv $SCRATCH/Archaea_v_Fungi/archaea_faa/out $SCRATCH/Archaea_v_Fungi/

cd $SCRATCH
rm arch_v_fung_script.tcl
