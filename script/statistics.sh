#! usr/bin/bash
PREFIX=Ath_root_NC

JOB=$(ls output/${PREFIX} | sort -r)
for J in $JOB;do
  echo "====>$J"

  pigz -dcf output/${PREFIX}/$J/${J}_align.sam.gz |
        grep -v "@" |
        tsv-filter --regex '6:^[0-9]+=$' |
        tsv-filter --str-eq 7:= |
        cut -f 1 | uniq > output/${PREFIX}/$J/${J}_match.tsv

  perl script/pre_remove_statistics.pl ../data/Ath_root_NC/R1.fq.gz output/${PREFIX}/$J/${J}_match.tsv \
  output/${PREFIX}/$J/reads_info.tsv \
  output/${PREFIX}/$J/remove_info.tsv
done
