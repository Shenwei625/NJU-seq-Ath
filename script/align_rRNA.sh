#! usr/bin/bash
for TISSUE in flower;do
    for PREFIX in Ath_${TISSUE}_NC Ath_${TISSUE}_1 Ath_${TISSUE}_2 Ath_${TISSUE}_3;do
        echo "===> align_$PREFIX"
        mkdir -p rRNA_conserve/output/${PREFIX}

        echo "===> $J conserve region align"
        mkdir -p rRNA_conserve/output/${PREFIX}/${J}

        bowtie2 -p 20 -a -t \
        --end-to-end -D 20 -R 3 \
        -N 0 -L 10 -i S,1,0.50 --np 0 \
        --xeq -x index/${J}/${J} \
        -1 ../data/${PREFIX}/R1.fq.gz -2 ../data/${PREFIX}/R2.fq.gz \
        -S rRNA_conserve/output/${PREFIX}/${J}/${J}_align.sam \
        2>&1 |
        tee rRNA_conserve/output/${PREFIX}/${J}/${J}.bowtie2.log 

        pigz -p 20 rRNA_conserve/output/${PREFIX}/${J}/${J}_align.sam
    done
done
