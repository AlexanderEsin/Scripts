	## Curl script for genome download ##

	ORG=organism
	cd ~/desktop/$ORG
	FILE=$ORG\_list.txt
	suffix=filetype



	while read line; do
		cd ~/desktop/$ORG/Bin;
		curl -O -m 10 --retry 10 --retry-delay 1 ftp://ftp.ncbi.nlm.nih.gov/genomes/all/$line/$line\_$suffix
	done < $FILE
	