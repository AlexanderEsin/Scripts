#!/bin/sh
#PBS -l walltime=50:00:00
#PBS -l mem=16gb
#PBS -l ncpus=64

module load blast+/2.2.28

####

org=FungArch_v_Bact
query_dir=fungarch_faa
db_dir=bacteria_db

####

cd $SCRATCH/$org/$query_dir
mkdir -p DB
mkdir -p out
mkdir -p timing

cd $SCRATCH/$org/$db_dir

dblist=?*.faa

for db in $dblist
do
	makeblastdb -in $db -dbtype prot -out $SCRATCH/$org/$query_dir/DB/$db.db
done


echo "
set org FungArch_v_Bact
set query_dir fungarch_faa
set db_dir bacteria_db

cd $SCRATCH/\$org/\$db_dir
set unsortlist [exec ls -s -S]
set header [string first \\n \$unsortlist]
set unsortlist [string range \$unsortlist \$header end]
set unsortlist [string trim \$unsortlist]

set pata {.+?[0-9]\s+?}
regsub -all -line \$pata \$unsortlist {} sortlist
puts \$sortlist
set dblist \$sortlist
cd $SCRATCH/\$org/\$query_dir
set glist [lrange [glob *.faa] [lindex \$argv 0] [lindex \$argv 1]]

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
	set out [open $SCRATCH/\$org/\$query_dir/timing/\$db.txt a]
	puts \$out \$err
	close \$out
}

puts DONE

exit
" > $SCRATCH/$org\_script.tcl

chmod +x $SCRATCH/$org\_script.tcl

tclsh8.5 $SCRATCH/$org\_script.tcl 0 5 &
tclsh8.5 $SCRATCH/$org\_script.tcl 6 11 &
tclsh8.5 $SCRATCH/$org\_script.tcl 12 17 &
tclsh8.5 $SCRATCH/$org\_script.tcl 18 23 &
tclsh8.5 $SCRATCH/$org\_script.tcl 24 29 &
tclsh8.5 $SCRATCH/$org\_script.tcl 30 35 &
tclsh8.5 $SCRATCH/$org\_script.tcl 36 41 &
tclsh8.5 $SCRATCH/$org\_script.tcl 42 47 &
tclsh8.5 $SCRATCH/$org\_script.tcl 48 53 &
tclsh8.5 $SCRATCH/$org\_script.tcl 54 59 &
tclsh8.5 $SCRATCH/$org\_script.tcl 60 65 &
tclsh8.5 $SCRATCH/$org\_script.tcl 66 71 &
tclsh8.5 $SCRATCH/$org\_script.tcl 72 77 &
tclsh8.5 $SCRATCH/$org\_script.tcl 78 83 &
tclsh8.5 $SCRATCH/$org\_script.tcl 84 89 &
tclsh8.5 $SCRATCH/$org\_script.tcl 90 95 &
tclsh8.5 $SCRATCH/$org\_script.tcl 96 101 &
tclsh8.5 $SCRATCH/$org\_script.tcl 102 107 &
tclsh8.5 $SCRATCH/$org\_script.tcl 108 113 &
tclsh8.5 $SCRATCH/$org\_script.tcl 114 119 &
tclsh8.5 $SCRATCH/$org\_script.tcl 120 125 &
tclsh8.5 $SCRATCH/$org\_script.tcl 126 131 &
tclsh8.5 $SCRATCH/$org\_script.tcl 132 137 &
tclsh8.5 $SCRATCH/$org\_script.tcl 138 143 &
tclsh8.5 $SCRATCH/$org\_script.tcl 144 149 &
tclsh8.5 $SCRATCH/$org\_script.tcl 150 155 &
tclsh8.5 $SCRATCH/$org\_script.tcl 156 161 &
tclsh8.5 $SCRATCH/$org\_script.tcl 162 167 &
tclsh8.5 $SCRATCH/$org\_script.tcl 168 173 &
tclsh8.5 $SCRATCH/$org\_script.tcl 174 179 &
tclsh8.5 $SCRATCH/$org\_script.tcl 180 185 &
tclsh8.5 $SCRATCH/$org\_script.tcl 186 191 &
tclsh8.5 $SCRATCH/$org\_script.tcl 192 197 &
tclsh8.5 $SCRATCH/$org\_script.tcl 198 203 &
tclsh8.5 $SCRATCH/$org\_script.tcl 204 209 &
tclsh8.5 $SCRATCH/$org\_script.tcl 210 215 &
tclsh8.5 $SCRATCH/$org\_script.tcl 216 221 &
tclsh8.5 $SCRATCH/$org\_script.tcl 222 227 &
tclsh8.5 $SCRATCH/$org\_script.tcl 228 233 &
tclsh8.5 $SCRATCH/$org\_script.tcl 234 239 &
tclsh8.5 $SCRATCH/$org\_script.tcl 240 245 &
tclsh8.5 $SCRATCH/$org\_script.tcl 246 251 &
tclsh8.5 $SCRATCH/$org\_script.tcl 252 257 &
tclsh8.5 $SCRATCH/$org\_script.tcl 258 263 &
tclsh8.5 $SCRATCH/$org\_script.tcl 264 269 &
tclsh8.5 $SCRATCH/$org\_script.tcl 270 275 &
tclsh8.5 $SCRATCH/$org\_script.tcl 276 281 &
tclsh8.5 $SCRATCH/$org\_script.tcl 282 287 &
tclsh8.5 $SCRATCH/$org\_script.tcl 288 293 &
tclsh8.5 $SCRATCH/$org\_script.tcl 294 299 &
tclsh8.5 $SCRATCH/$org\_script.tcl 300 305 &
tclsh8.5 $SCRATCH/$org\_script.tcl 306 311 &
tclsh8.5 $SCRATCH/$org\_script.tcl 312 317 &
tclsh8.5 $SCRATCH/$org\_script.tcl 318 323 &
tclsh8.5 $SCRATCH/$org\_script.tcl 324 329 &
tclsh8.5 $SCRATCH/$org\_script.tcl 330 335 &
tclsh8.5 $SCRATCH/$org\_script.tcl 336 341 &
tclsh8.5 $SCRATCH/$org\_script.tcl 342 347 &
tclsh8.5 $SCRATCH/$org\_script.tcl 348 353 &
tclsh8.5 $SCRATCH/$org\_script.tcl 354 359 &
tclsh8.5 $SCRATCH/$org\_script.tcl 360 365 &
tclsh8.5 $SCRATCH/$org\_script.tcl 366 372 &
tclsh8.5 $SCRATCH/$org\_script.tcl 373 379 &
tclsh8.5 $SCRATCH/$org\_script.tcl 380 386 &


wait

rm -r $SCRATCH/$org/$query_dir/DB
mv $SCRATCH/$org/$query_dir/out $SCRATCH/$org/

cd $SCRATCH
rm $org\_script.tcl
