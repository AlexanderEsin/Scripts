#!/bin/sh
#PBS -l walltime=50:00:00
#PBS -l mem=16gb
#PBS -l ncpus=64

module load blast+/2.2.28

cd $SCRATCH/Fungi_v_Archaea/fungi_faa
mkdir -p DB
mkdir -p out
mkdir -p timing

cd $SCRATCH/Fungi_v_Archaea/archaea_db

dblist=?*.faa

for db in $dblist
do
	makeblastdb -in $db -dbtype prot -out $SCRATCH/Fungi_v_Archaea/fungi_faa/DB/$db.db
done


echo "cd $SCRATCH/Fungi_v_Archaea/archaea_db
set unsortlist [exec ls -s -S]
set header [string first \\n \$unsortlist]
set unsortlist [string range \$unsortlist \$header end]
set unsortlist [string trim \$unsortlist]

set pata {.+?[0-9]\s+?}
regsub -all -line \$pata \$unsortlist {} sortlist
puts \$sortlist
set dblist [lrange \$sortlist [lindex \$argv 0] [lindex \$argv 1]]
cd $SCRATCH/Fungi_v_Archaea/fungi_faa
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
			catch {exec /apps/ncbi-blast/2.2.28/bin/blastp -query \$g -db DB/\$db.db -out out/\$x\&\$y.tsv -evalue 1e-10 -outfmt 7 -max_target_seqs 1 -max_hsps_per_subject 1 -num_threads 6}
			set systime [clock seconds]
			append err "\$g Finish: [clock format \$systime -format %H:%M:%S]\\n\\n"
		}
	}
	set out [open $SCRATCH/Fungi_v_Archaea/fungi_faa/timing/\$db.txt a]
	puts \$out \$err
	close \$out
}" > $SCRATCH/fung_v_arch_script.tcl

chmod +x $SCRATCH/fung_v_arch_script.tcl

tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 0 2 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 3 5 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 6 8 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 9 12 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 13 16 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 17 20 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 21 24 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 25 28 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 29 32 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 33 36 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 37 40 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 41 44 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 45 48 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 49 52 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 53 56 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 57 60 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 61 64 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 65 68 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 69 72 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 73 76 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 77 80 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 81 84 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 85 88 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 89 92 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 93 96 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 97 100 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 101 104 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 105 108 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 109 112 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 113 116 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 117 120 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 121 124 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 125 128 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 129 132 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 133 136 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 137 140 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 141 144 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 145 148 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 149 152 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 153 156 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 157 160 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 161 164 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 165 168 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 169 172 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 173 176 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 177 180 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 181 184 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 185 188 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 189 192 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 193 196 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 197 200 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 201 204 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 205 208 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 209 212 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 213 216 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 217 220 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 221 224 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 225 228 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 229 232 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 233 236 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 237 240 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 241 244 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 245 248 &
tclsh8.5 $SCRATCH/fung_v_arch_script.tcl 249 252 &

wait

rm -r $SCRATCH/Fungi_v_Archaea/fungi_faa/DB
mv $SCRATCH/Fungi_v_Archaea/fungi_faa/out $SCRATCH/Fungi_v_Archaea/

cd $SCRATCH
rm fung_v_arch_script.tcl
