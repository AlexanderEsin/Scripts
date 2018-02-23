# In one session start the master
/usr/local/bin/hmmpgmd --master --cport 51500 --wport 51501 --hmmdb /Users/aesin/Desktop/emapper/data/hmmdb_levels/bact_50/bact_50.hmm

# In the other start the worker
/usr/local/bin/hmmpgmd --worker localhost --wport 51501 --cpu 22