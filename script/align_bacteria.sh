#! usr/bin/bash
# bash align.sh PREFIX
PREFIX=$1

mkdir -p rRNA_tRNA/output/${PREFIX}

find index/rRNA_tRNA -maxdepth 1 -mindepth 1 -type d |
    cut -d "/" -f 3 |
      parallel -j 20 --line-buffer "
        echo '===>$PREFIX {1}'
        mkdir -p rRNA_tRNA/output/${PREFIX}/{1}

        bowtie2 -a -t \
          --end-to-end -D 20 -R 3 \
          -N 0 -L 10 -i S,1,0.50 --np 0 \
          --xeq -x index/rRNA_tRNA/{1}/{1} \
          -1 rRNA_conserve/filter/${PREFIX}/R1_conserve.fq.gz -2 rRNA_conserve/filter/${PREFIX}/R2_conserve.fq.gz \
          -S rRNA_tRNA/output/${PREFIX}/{1}/{1}_align.sam \
          2>&1 |
        tee rRNA_tRNA/output/${PREFIX}/{1}/{1}.bowtie2.log

        pigz rRNA_tRNA/output/${PREFIX}/{1}/{1}_align.sam
      "
