#! usr/bin/bash
# bash align_mRNA.sh PREFIX JOB
PREFIX=$1
JOB=$2

cat $2 |
  parallel -j 20 --line-buffer "
    echo '===>$PREFIX {1}'
    mkdir -p mRNA/output/${PREFIX}/{1}

    bowtie2 -a -t \
    --end-to-end -D 20 -R 3 \
    -N 0 -L 10 -i S,1,0.50 --np 0 \
    --xeq -x index/mRNA/{1}/{1} \
    -1 rRNA_tRNA/filter/${PREFIX}/R1_conserve_rRNA.fq.gz -2 rRNA_tRNA/filter/${PREFIX}/R2_conserve_rRNA.fq.gz \
    -S mRNA/output/${PREFIX}/{1}/{1}_align.sam \
    2>&1 |
    tee mRNA/output/${PREFIX}/{1}/{1}.bowtie2.log

    pigz mRNA/output/${PREFIX}/{1}/{1}_align.sam
  "
