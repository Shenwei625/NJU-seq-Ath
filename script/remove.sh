#! usr/bin/bash
for PREFIX in Ath_flower_1 Ath_flower_2 Ath_flower_3;do
echo "===> ${PREFIX}"
mkdir -p rRNA_conserve/filter/${PREFIX}

pigz -dcf ../data/${PREFIX}/R1.fq.gz | grep "@" |
    cut -d " " -f 1 | sed 's/^@//g' |
    grep -v -w -f <(sed '1d' rRNA_conserve/output/${PREFIX}/total_remove_reads_info.tsv | cut -f 1 ) > rRNA_conserve/filter/${PREFIX}/keep.lst

for J in R1 R2;do
  echo "====> $J"
  seqkit grep -j 4 -f rRNA_conserve/filter/${PREFIX}/keep.lst ../data/${PREFIX}/${J}.fq.gz > rRNA_conserve/filter/${PREFIX}/${J}_conserve.fq
done

done

