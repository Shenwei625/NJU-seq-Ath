#! usr/bin/bash

PREFIX=$1

gzip -dcf stage_align/${PREFIX}/mrna.new.sam.gz |
  parallel --pipe --block 100M --no-run-if-empty --linebuffer --keep-order -j 24 '
    awk '\''$6!="*"&&$7=="="{print $1 "\t" $3 "\t" $4 "\t" $6 "\t" $10}
    '\'' |
    perl NJU_seq/mrna_analysis/multimatch_judge.pl
  ' | perl NJU_seq/mrna_analysis/multimatch_judge.pl \
  >temp/${PREFIX}/mrna.out.tmp
