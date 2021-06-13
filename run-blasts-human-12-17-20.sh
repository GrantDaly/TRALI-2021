#!/bin/sh
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load blast+/2.6.0



run_blast () {
if [[ -p /scratch/blast-fifo ]]; then
rm /scratch/blast-fifo
fi

#make out pipe
if [[ -p /scratch/out-pipe ]]; then
rm /scratch/out-pipe
fi

mkfifo /scratch/out-pipe
mkfifo /scratch/blast-fifo

echo "unzipping input"
zcat $in_dir/$in_fasta | sed -n '1~4s/^@/>/p;2~4p' > /scratch/blast-fifo &

#set up gzip
gzip -9 -c < /scratch/out-pipe > $OUTNAME.gz &
echo "starting blast"
date
blastn -db "Homo_Sapiens" \
-query /scratch/blast-fifo \
-word_size 9 \
-reward 1 \
-penalty -1 \
-gapopen 2 \
-gapextend 2 \
-out /scratch/out-pipe \
-num_threads $THREADS \
-outfmt 15
#-outfmt "6 qseqid qstart qend sseqid sstart send sstrand pident length mismatch gapopend evalue bitscore stitle qseq sseq"
echo "done with blast"
date
}

THREADS=16
in_dir=/home/usagtd/mitochondrial-NUMTs
in_fasta=mitochondrial.fasta.gz
OUTNAME=rCStoGRCh38-12-17-20.json
run_blast
