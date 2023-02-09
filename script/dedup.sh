#! usr/bin/bash

PREFIX=$1

cat temp/${PREFIX}/mrna.out.tmp |
  parallel --pipe --block 100M --no-run-if-empty --linebuffer --keep-order -j 24 '
    perl NJU_seq/mrna_analysis/dedup.pl \
      --refstr "Parent=transcript:" \
      --transid "AT" \
      --info data/ath_exon.info
  ' |
  perl NJU_seq/mrna_analysis/dedup.pl \
    --refstr "Parent=transcript:" \
    --transid "AT" \
    --info data/ath_exon.info \
    >temp/${PREFIX}/mrna.dedup.tmp
