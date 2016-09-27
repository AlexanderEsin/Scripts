#!/bin/sh
#PBS -l walltime=50:00:00
#PBS -l mem=16gb
#PBS -l ncpus=64

module load blast+/2.2.28

cd $SCRATCH/Archaea_v_Fungi/archaea_faa
mkdir -p DB
mkdir -p out2
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
			catch {exec /apps/ncbi-blast/2.2.28/bin/blastp -query \$g -db DB/\$db.db -out out2/\$x\&\$y.tsv -evalue 1e-10 -outfmt 7 -max_target_seqs 1 -max_hsps_per_subject 1 -num_threads 6}
			set systime [clock seconds]
			append err "\$g Finish: [clock format \$systime -format %H:%M:%S]\\n\\n"
		}
	}
	set out [open $SCRATCH/Archaea_v_Fungi/archaea_faa/timing/\$db.txt a]
	puts \$out \$err
	close \$out
}" > $SCRATCH/arch_v_fung_script.tcl

chmod +x $SCRATCH/arch_v_fung_script.tcl

tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 0 1 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 2 3 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 4 5 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 6 7 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 8 9 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 10 11 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 12 13 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 14 15 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 16 17 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 18 19 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 20 21 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 22 23 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 24 25 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 26 27 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 28 29 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 30 31 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 32 33 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 34 35 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 36 37 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 38 39 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 40 41 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 42 43 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 44 45 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 46 47 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 48 49 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 50 51 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 52 53 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 54 55 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 56 57 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 58 59 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 60 61 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 62 63 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 64 65 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 66 67 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 68 69 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 70 71 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 72 73 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 74 75 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 76 77 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 78 79 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 80 81 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 82 83 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 84 85 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 86 87 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 88 89 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 90 91 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 92 93 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 94 95 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 96 97 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 98 99 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 100 101 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 102 103 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 104 105 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 106 107 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 108 109 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 110 111 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 112 113 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 114 115 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 116 118 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 119 121 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 122 124 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 125 127 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 128 130 &
tclsh8.5 $SCRATCH/arch_v_fung_script.tcl 131 133 &

wait

rm -r $SCRATCH/Archaea_v_Fungi/archaea_faa/DB
mv $SCRATCH/Archaea_v_Fungi/archaea_faa/out2 $SCRATCH/Archaea_v_Fungi/

cd $SCRATCH
rm arch_v_fung_script.tcl
