while IFS='' read -r line || [[ -n "$line" ]]; do
    fastq-dump --outdir /Users/aesin/Desktop/operonExpr/SRA_files $line
done < "$1"