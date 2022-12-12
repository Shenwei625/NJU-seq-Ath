#! usr/bin/bash
# bash statistics.sh dir PREFIX

PREFIX=$2

ls $1/output/${PREFIX} |
  parallel -j 20 --linebuffer "
  echo '====> {1}'

  pigz -dcf $1/output/${PREFIX}/{1}/{1}_align.sam.gz |
        grep -v "@" |
        tsv-filter --regex '6:^[0-9]+=$' |
        tsv-filter --str-eq 7:= |
        cut -f 1,10 | tsv-uniq -f 1 > $1/output/${PREFIX}/{1}/{1}_match.tsv

  perl script/pre_remove_statistics.pl ../data/${PREFIX}/R1.fq.gz $1/output/${PREFIX}/{1}/{1}_match.tsv \
  $1/output/${PREFIX}/{1}/reads_info.tsv \
  $1/output/${PREFIX}/{1}/remove_info.tsv
"
