import os
import glob
import re
import sys

f = open(os.devnull, 'w')

markdir = "/users/aesin/desktop/test_blast/input/marked_faa"

os.chdir(markdir)
dblist = glob.glob('*faa')

if not os.path.exists('DB'):
    os.makedirs('DB')

if not os.path.exists('Out'):
    os.makedirs('Out')

for db in dblist:
	os.system('makeblastdb -in %s -dbtype prot -out DB/%s.db'%(db,db))
	glist = glob.glob('*faa')
	for g in glist:
		if g == db:
			print 'Skip'
		else:
			x = re.sub('_protein.faa', "", g)
			y = re.sub('_protein.faa', "", db)
			print('%s vs %s' %(g,db))
			os.system('blastp -query %s -db DB/%s.db -out Out/%s\&%s.tsv -evalue 1e-10 -outfmt 6 -max_target_seqs 1 -max_hsps 1 -num_threads 3' % (g,db,x,y))