#! usr/bin/bash
PREFIX=Ath_flower_NC

JOB=$(ls output/${PREFIX})

echo -e "name\tTotal_reads_number\tMatch_reads_number\tMatch_proportion\tAverage_CG_content\tAverage_reads_Length" > output/${PREFIX}/align_statistics.tsv

for J in $JOB;do
  echo "====> $J"

  cat output/${PREFIX}/$J/remove_info.tsv |  
    datamash transpose > output/${PREFIX}/$J/tem&&
    mv output/${PREFIX}/$J/tem output/${PREFIX}/$J/remove_info.tsv

  # strain name
  (echo -e "bacteria_name\t$J" && cat output/${PREFIX}/$J/remove_info.tsv) > output/${PREFIX}/$J/tem &&
    mv output/${PREFIX}/$J/tem output/${PREFIX}/$J/remove_info.tsv

  # Average CG content
  CG_AVERAGE=$(cat output/${PREFIX}/$J/reads_info.tsv | perl -e'
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
  echo -e "Average_CG_content\t$CG_AVERAGE" >> output/${PREFIX}/$J/remove_info.tsv

  # Average reads length
  length_AVERAGE=$(cat output/${PREFIX}/$J/reads_info.tsv | perl -e'
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
  echo -e "Average_reads_length\t$length_AVERAGE" >> output/${PREFIX}/$J/remove_info.tsv

  cat output/${PREFIX}/$J/remove_info.tsv |  
    datamash transpose > output/${PREFIX}/$J/tem&&
    mv output/${PREFIX}/$J/tem output/${PREFIX}/$J/remove_info.tsv

  # Length distribution
  (echo -e "Reads_length\tNumber" && cut -f 3 output/${PREFIX}/$J/reads_info.tsv | sort -n |
    uniq --count | perl -ne'
      s/^\s+//;
      print "$_";
    ' | tr " " "\t" |
    tsv-select --fields 2,1) > output/${PREFIX}/$J/reads_length_distribution.tsv

  # merge
  sed '1d' output/${PREFIX}/$J/remove_info.tsv >> output/${PREFIX}/align_statistics.tsv
done
