#! usr/bin/bash
PREFIX=Ath_flower_NC

ls ./output/${PREFIX} |
  parallel -j 24 --linebuffer '
  echo "====> {1}"

  pigz -dcf ./output/Ath_flower_NC/{1}/{1}_align.sam.gz |
        grep -v "@" |
        tsv-filter --regex '6:^[0-9]+=$' |
        tsv-filter --str-eq 7:= |
        cut -f 1,10 | uniq > ./output/Ath_flower_NC/{1}/{1}_match.tsv

  perl script/pre_remove_statistics.pl ../data/Ath_flower_NC/R1.fq.gz output/Ath_flower_NC/{1}/{1}_match.tsv \
  output/Ath_flower_NC/{1}/reads_info.tsv \
  output/Ath_flower_NC/{1}/remove_info.tsv
'