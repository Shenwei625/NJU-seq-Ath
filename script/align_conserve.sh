#! usr/bin/bash
# bash align_conserve.sh PREFIX
PREFIX=$1

for J in 16S 23S;do
    echo "===> $J conserve region align"
    mkdir -p rRNA_conserve/output/${PREFIX}/${J}

    bowtie2 -p 24 -a -t \
    --end-to-end -D 20 -R 3 \
    -N 0 -L 10 -i S,1,0.50 --np 0 \
    --xeq -x index/${J}/${J} \
    -1 ../data/${PREFIX}/R1.fq.gz -2 ../data/${PREFIX}/R2.fq.gz \
    -S rRNA_conserve/output/${PREFIX}/${J}/${J}_align.sam \
    2>&1 |
    tee rRNA_conserve/output/${PREFIX}/${J}/${J}.bowtie2.log 

    pigz -p 24 rRNA_conserve/output/${PREFIX}/${J}/${J}_align.sam
done
