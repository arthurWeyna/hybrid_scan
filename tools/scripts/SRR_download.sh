##This script downloads fastq files from genbank using fasterq-dump
#$1 is a file with ids to download in first column (such as uce_list)
#$2 is the directory where to dump fastq.gz files

ids=`cut -f1 -d" " $1 | grep -v "#"`
#ids="SRR5760503"
out=$2
threads=5

for id in $ids
do
	echo $id	
	if [ `ls $out | grep $id | wc -l` -lt 1 ];
	then
		fasterq-dump --outdir $out --temp "$out"/"$id"_sratmp --mem 5G --split-3 --threads $threads --skip-technical  --print-read-nr $id
		rm -rf "$out"/"$id"_sratmp
		pigz -p"$threads" "$out"/"$id"*.fastq
	fi
done
