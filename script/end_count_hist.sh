#! usr/bin/bash
# bash end_count_hist.sh ${TISSUE} rrna_5-8s/18s/25s/mrna

TISSUE=$1
TARGET=$2

# NC
echo -e "Site\tbase\tEnd_count" > visuliaztion/End_count_hist/${TISSUE}_${TARGET}_NC.tsv
sed '$d' output/Ath_${TISSUE}_NC/${TARGET}.tsv | cut -f 1,2,4 >>  visuliaztion/End_count_hist/${TISSUE}_${TARGET}_NC.tsv
# TR
sed '$d' output/Ath_${TISSUE}_1/${TARGET}.tsv | cut -f 1,2,4 >>  visuliaztion/End_count_hist/${TISSUE}_${TARGET}_TR.tsv
tsv-join --filter-file output/Ath_${TISSUE}_2/${TARGET}.tsv --key-fields 1 --append-fields 4 visuliaztion/End_count_hist/${TISSUE}_${TARGET}_TR.tsv > tem&&
    mv tem visuliaztion/End_count_hist/${TISSUE}_${TARGET}_TR.tsv
tsv-join --filter-file output/Ath_${TISSUE}_3/${TARGET}.tsv --key-fields 1 --append-fields 4 visuliaztion/End_count_hist/${TISSUE}_${TARGET}_TR.tsv > tem&&
    mv tem visuliaztion/End_count_hist/${TISSUE}_${TARGET}_TR.tsv

cat visuliaztion/End_count_hist/${TISSUE}_${TARGET}_TR.tsv | perl -a -F"\t" -ne '
  chomp;
  my $END = ( ( @F[2] + @F[3] + @F[4] ) / 3 );
  print "$_\t";
  printf "%.0f\n", $END;
' > tem&&
mv tem visuliaztion/End_count_hist/${TISSUE}_${TARGET}_TR.tsv

(echo -e "Site\tbase\tTR1_end_count\tTR2_end_count\tTR3_end_count\tEnd_count" && cat visuliaztion/End_count_hist/${TISSUE}_${TARGET}_TR.tsv) > tem&&
    mv tem visuliaztion/End_count_hist/${TISSUE}_${TARGET}_TR.tsv

Rscript remove_bacteria/script/end_count_hist.r ${TISSUE} ${TARGET}
rm visuliaztion/End_count_hist/${TISSUE}_${TARGET}_TR.tsv
rm visuliaztion/End_count_hist/${TISSUE}_${TARGET}_NC.tsv
