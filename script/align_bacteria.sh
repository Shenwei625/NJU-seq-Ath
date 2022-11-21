#! usr/bin/bash
for TISSUE in flower stem leaf;do
  for PREFIX in Ath_${TISSUE}_NC Ath_${TISSUE}_1 Ath_${TISSUE}_2 Ath_${TISSUE}_3;do
    THREAD=24
      find ASSEMBLY -maxdepth 1 -mindepth 1 -type d |
      cut -d "/" -f 2 |
      parallel -j 12 --line-buffer "
        echo '===>{1}'
        mkdir -p output/${PREFIX}/{1}

        bowtie2 -p ${THREAD} -a -t \
          --end-to-end -D 20 -R 3 \
          -N 0 -L 10 -i S,1,0.50 --np 0 \
          --xeq -x index/{1}/{1} \
          -1 ../data/${PREFIX}/R1.fq.gz -2 ../data/${PREFIX}/R2.fq.gz \
          -S output/${PREFIX}/{1}/{1}_align.sam \
          2>&1 |
        tee output/${PREFIX}/{1}/{1}.bowtie2.log

        pigz -p ${THREAD} output/${PREFIX}/{1}/{1}_align.sam
      "
  done
done