#!/bin/sh
#PBS -l walltime=50:00:00
#PBS -l mem=12gb
#PBS -l ncpus=32

module load blast+/2.2.28

####

org=Bact_v_FungArch
query_dir=bacteria_faa
db_dir=fungarch_db

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
set org Bact_v_FungArch
set query_dir bacteria_faa
set db_dir fungarch_db

cd $SCRATCH/\$org/\$db_dir
set unsortlist [exec ls -s -S]
set header [string first \\n \$unsortlist]
set unsortlist [string range \$unsortlist \$header end]
set unsortlist [string trim \$unsortlist]

set pata {.+?[0-9]\s+?}
regsub -all -line \$pata \$unsortlist {} sortlist
puts \$sortlist
set dblist [lrange \$sortlist [lindex \$argv 0] [lindex \$argv 1]]
cd $SCRATCH/\$org/\$query_dir
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
	set out [open $SCRATCH/\$org/\$query_dir/timing/\$db.txt a]
	puts \$out \$err
	close \$out
}

puts DONE

exit
" > $SCRATCH/$org\_script.tcl

chmod +x $SCRATCH/$org\_script.tcl

tclsh8.5 $SCRATCH/$org\_script.tcl 0 11 &
tclsh8.5 $SCRATCH/$org\_script.tcl 12 23 &
tclsh8.5 $SCRATCH/$org\_script.tcl 24 35 &
tclsh8.5 $SCRATCH/$org\_script.tcl 36 47 &
tclsh8.5 $SCRATCH/$org\_script.tcl 48 59 &
tclsh8.5 $SCRATCH/$org\_script.tcl 60 71 &
tclsh8.5 $SCRATCH/$org\_script.tcl 72 83 &
tclsh8.5 $SCRATCH/$org\_script.tcl 84 95 &
tclsh8.5 $SCRATCH/$org\_script.tcl 96 107 &
tclsh8.5 $SCRATCH/$org\_script.tcl 108 119 &
tclsh8.5 $SCRATCH/$org\_script.tcl 120 131 &
tclsh8.5 $SCRATCH/$org\_script.tcl 132 143 &
tclsh8.5 $SCRATCH/$org\_script.tcl 144 155 &
tclsh8.5 $SCRATCH/$org\_script.tcl 156 167 &
tclsh8.5 $SCRATCH/$org\_script.tcl 168 179 &
tclsh8.5 $SCRATCH/$org\_script.tcl 180 191 &
tclsh8.5 $SCRATCH/$org\_script.tcl 192 203 &
tclsh8.5 $SCRATCH/$org\_script.tcl 204 215 &
tclsh8.5 $SCRATCH/$org\_script.tcl 216 227 &
tclsh8.5 $SCRATCH/$org\_script.tcl 228 239 &
tclsh8.5 $SCRATCH/$org\_script.tcl 240 251 &
tclsh8.5 $SCRATCH/$org\_script.tcl 252 263 &
tclsh8.5 $SCRATCH/$org\_script.tcl 264 275 &
tclsh8.5 $SCRATCH/$org\_script.tcl 276 287 &
tclsh8.5 $SCRATCH/$org\_script.tcl 288 299 &
tclsh8.5 $SCRATCH/$org\_script.tcl 300 311 &
tclsh8.5 $SCRATCH/$org\_script.tcl 312 323 &
tclsh8.5 $SCRATCH/$org\_script.tcl 324 335 &
tclsh8.5 $SCRATCH/$org\_script.tcl 336 347 &
tclsh8.5 $SCRATCH/$org\_script.tcl 348 360 &
tclsh8.5 $SCRATCH/$org\_script.tcl 361 373 &
tclsh8.5 $SCRATCH/$org\_script.tcl 374 386 &



wait

rm -r $SCRATCH/$org/$query_dir/DB
mv $SCRATCH/$org/$query_dir/out $SCRATCH/$org/

cd $SCRATCH
rm $org\_script.tcl
