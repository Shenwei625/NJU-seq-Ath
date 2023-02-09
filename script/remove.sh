#! usr/bin/bash
# bash remove.sh dir PREFIX

dr=$1
PREFIX=$2

mkdir -p ${dr}/filter/${PREFIX}

pigz -dcf ../data/${PREFIX}/R1.fq.gz | grep "@" |
    cut -d " " -f 1 | sed 's/^@//g' |
    grep -v -w -f <(sed '1d' ${dr}/output/${PREFIX}/total_remove_reads_info.tsv | cut -f 1 ) > ${dr}/filter/${PREFIX}/keep.lst

for J in R1 R2;do
  echo "====> $J"
  seqkit grep -j 4 -f ${dr}/filter/${PREFIX}/keep.lst ../data/${PREFIX}/${J}.fq.gz > ${dr}/filter/${PREFIX}/${J}_conserve.fq

  pigz -p 4 ${dr}/filter/${PREFIX}/${J}_conserve.fq
done
