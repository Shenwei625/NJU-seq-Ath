#! usr/bin/bash
# bash stage_align.sh stage/remove DIR PREFIX rrna/protein_coding
DIR=$2
PREFIX=$3

mkdir -p ${1}_align/${PREFIX}/${DIR}

if [ $1 == stage ]; then
    DIR2=filter
else
    DIR2=remove
fi
echo $DIR2

if [ $4 == rrna ]; then
    RNA=rrna
else
    RNA=mrna
fi
echo $RNA

bowtie2 -p 24 -a -t \
    --end-to-end -D 20 -R 3 \
    -N 0 -L 10 -i S,1,0.50 --np 0 \
    --xeq -x index/ath_${4} \
    -1 remove_bacteria/${DIR}/${DIR2}/${PREFIX}/R1*.fq.gz -2 remove_bacteria/${DIR}/${DIR2}/${PREFIX}/R2*.fq.gz \
    -S ${1}_align/${PREFIX}/${DIR}/${RNA}_align.sam \
    2>&1 |
    tee ${1}_align/${PREFIX}/${DIR}/${RNA}.bowtie2.log

pigz -p 24 ${1}_align/${PREFIX}/${DIR}/${RNA}_align.sam
