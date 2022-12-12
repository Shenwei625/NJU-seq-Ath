#! usr/bin/bash
# bash merge_plot.sh dir PREFIX
PREFIX=$2

JOB=$(ls $1/output/${PREFIX})

echo -e "name\tTotal_reads_number\tMatch_reads_number\tMatch_proportion\tAverage_CG_content\tAverage_reads_Length" > $1/output/${PREFIX}/align_statistics.tsv

for J in $JOB;do
  echo "====> calculate $J"

  cat $1/output/${PREFIX}/$J/remove_info.tsv |  
    datamash transpose > $1/output/${PREFIX}/$J/tem&&
    mv $1/output/${PREFIX}/$J/tem $1/output/${PREFIX}/$J/remove_info.tsv

  # strain name
  (echo -e "bacteria_name\t$J" && cat $1/output/${PREFIX}/$J/remove_info.tsv) > $1/output/${PREFIX}/$J/tem &&
    mv $1/output/${PREFIX}/$J/tem $1/output/${PREFIX}/$J/remove_info.tsv

  # Average CG content
  CG_AVERAGE=$(cat $1/output/${PREFIX}/$J/reads_info.tsv | perl -e'
              my $total_CG;
              my $n;
              while (<>) {
                if (/\S+\t\S+\t\S+\t(\S+)/) {
                  $total_CG = $total_CG + $1;
                  $n++;
                }
              }
              my $CG_average = ( ( $total_CG / $n ) * 100 );
              printf "%.3f%%", $CG_average;
              ')
  echo -e "Average_CG_content\t$CG_AVERAGE" >> $1/output/${PREFIX}/$J/remove_info.tsv

  # Average reads length
  length_AVERAGE=$(cat $1/output/${PREFIX}/$J/reads_info.tsv | perl -e'
                my $total_length;
                my $n;
                while (<>) {
                  if (/\S+\t\S+\t(\S+)\t\S+/) {
                    $total_length = $total_length + $1;
                    $n++;
                  }
                }
                my $length_average = ( $total_length / $n );
                printf "%.2f", $length_average;
                ')
  echo -e "Average_reads_length\t$length_AVERAGE" >> $1/output/${PREFIX}/$J/remove_info.tsv

  cat $1/output/${PREFIX}/$J/remove_info.tsv |  
    datamash transpose > $1/output/${PREFIX}/$J/tem&&
    mv $1/output/${PREFIX}/$J/tem $1/output/${PREFIX}/$J/remove_info.tsv

  # Length distribution
  (echo -e "Reads_length\tNumber" && cut -f 3 $1/output/${PREFIX}/$J/reads_info.tsv | sed '1d' | sort -n |
    uniq --count | perl -ne'
      s/^\s+//;
      print "$_";
    ' | tr " " "\t" |
    tsv-select --fields 2,1) > $1/output/${PREFIX}/$J/reads_length_distribution.tsv

  sed '1d' $1/output/${PREFIX}/$J/remove_info.tsv >> $1/output/${PREFIX}/align_statistics.tsv
done

cat $1/output/${PREFIX}/align_statistics.tsv | keep-header -- sort -nr -k4,4 > $1/output/${PREFIX}/align_tem&&
  mv $1/output/${PREFIX}/align_tem $1/output/${PREFIX}/align_statistics.tsv

# extract
for J in $JOB;do
  echo "===> merge $J"
  cat $1/output/${PREFIX}/$J/reads_info.tsv >> $1/output/${PREFIX}/total_remove_reads_info.tsv
done

tsv-uniq -f 1 $1/output/${PREFIX}/total_remove_reads_info.tsv > $1/output/${PREFIX}/total_tem&&
  mv $1/output/${PREFIX}/total_tem $1/output/${PREFIX}/total_remove_reads_info.tsv
