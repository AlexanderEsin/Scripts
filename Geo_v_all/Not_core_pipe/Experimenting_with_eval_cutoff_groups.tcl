time {
	###
	set direct /users/aesin/desktop/Geo_v_all/Group_fastas
	###

	proc openfile {fl} {
		global data
		set in [open $fl r]
		set data [read $in]
		close $in
		return
	}

	###########################################################################

	cd $direct

	file mkdir Eval_cutoff
	set pata {.+?\n\n+?}
	set i 1
	set groups [glob *.faa]
	set groups [lsort -dictionary $groups]
	set groups [lrange $groups 0 4]
	foreach group $groups {
		regsub {.faa} $group {} number
		file mkdir Eval_cutoff/$number

		openfile $group
		regsub -all {>} $data "\n>" glist
		set glist "$glist\n"
		regexp -line {>+?.+?(kaustophilus).+?\n+?} $glist comment_line
		regexp {>+?.+?\s+?} $comment_line id
		regexp $id$pata $glist gene

		set out [open $direct/Eval_cutoff/$number/query.faa w]
		puts $out [string trim $gene]
		close $out

		set out [open $direct/Eval_cutoff/$number/database.faa w]
		puts $out [string trim $data]
		close $out

		catch {exec makeblastdb -in $direct/Eval_cutoff/$number/database.faa -dbtype prot -out $direct/Eval_cutoff/$number/$number\_db.db}

		catch {exec blastp -query $direct/Eval_cutoff/$number/query.faa -db $direct/Eval_cutoff/$number/$number\_db.db -out $direct/Eval_cutoff/$number/$number\_output.tsv -evalue 1e-50 -outfmt "6 pident slen qcovs evalue" -max_target_seqs 10000 -max_hsps 1 -num_threads 8}

		puts "$i / [llength $groups]"
		incr i
	}
}
