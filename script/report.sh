#! usr/bin/bash
# bash report.sh PREFIX
PREFIX=$1

echo -e "Stage\treads\tproportion\trRNA_align_rate\tmRNA_align_rate" > visuliaztion/${PREFIX}_align_rate.tsv
# raw
export RAW=$(cat raw_output/${PREFIX}/rrna.bowtie2.log | 
  perl -ne '
    if (/^(\d+) reads; of these:/) {
      print "$1";
    }'
    )

cat raw_output/${PREFIX}/rrna.bowtie2.log | perl -ne '
  if (/^(\d+) reads; of these:/) {
    print "raw\t$1\t";
    my $proportion = ( ( $1 / $ENV{RAW} ) * 100 );
    printf "%.2f%%\t", $proportion;
  }
  if (/^(\S+) overall alignment rate/) {
    print "$1\t";
  }
' >> visuliaztion/${PREFIX}_align_rate.tsv

cat raw_output/${PREFIX}/mrna.bowtie2.log | perl -ne '
  if (/^(\S+) overall alignment rate/) {
    print "$1\n";
  }
' >> visuliaztion/${PREFIX}_align_rate.tsv

# different_stage
for DIR in rRNA_conserve rRNA_tRNA mRNA;do
  export i=$DIR
  cat stage_align/${PREFIX}/${DIR}/rrna.bowtie2.log | perl -ne '
    if (/^(\d+) reads; of these:/) {
      print "$ENV{i}_keep\t$1\t";
      my $proportion = ( ( $1 / $ENV{RAW} ) * 100 );
      printf "%.2f%%\t", $proportion;
    }
    if (/^(\S+) overall alignment rate/) {
      print "$1\t";
    } 
  ' >> visuliaztion/${PREFIX}_align_rate.tsv

  cat stage_align/${PREFIX}/${DIR}/mrna.bowtie2.log | perl -ne '
    if (/^(\S+) overall alignment rate/) {
      print "$1\n";
    } 
  ' >> visuliaztion/${PREFIX}_align_rate.tsv

  cat remove_align/${PREFIX}/${DIR}/rrna.bowtie2.log | perl -ne '
    if (/^(\d+) reads; of these:/) {
      print "$ENV{i}_remove\t$1\t";
      my $proportion = ( ( $1 / $ENV{RAW} ) * 100 );
      printf "%.2f%%\t", $proportion;
    }
    if (/^(\S+) overall alignment rate/) {
      print "$1\t";
    } 
  ' >> visuliaztion/${PREFIX}_align_rate.tsv

    cat remove_align/${PREFIX}/${DIR}/mrna.bowtie2.log | perl -ne '
    if (/^(\S+) overall alignment rate/) {
      print "$1\n";
    } 
  ' >> visuliaztion/${PREFIX}_align_rate.tsv
done

echo "<div align=center>" > visuliaztion/${PREFIX}_align_info.md
echo -e "# ${PREFIX}" >> visuliaztion/${PREFIX}_align_info.md

# raw
echo -e "## raw\n" >> visuliaztion/${PREFIX}_align_info.md
head -n 2 visuliaztion/${PREFIX}_align_rate.tsv | 
    mlr --itsv --omd cat >> visuliaztion/${PREFIX}_align_info.md
echo -e "\n" >> visuliaztion/${PREFIX}_align_info.md

# rRNA_conserve
for DIR in rRNA_conserve;do
    echo -e "## ${DIR}\n" >> visuliaztion/${PREFIX}_align_info.md
    (head -n 1 visuliaztion/${PREFIX}_align_rate.tsv && grep "^${DIR}" visuliaztion/${PREFIX}_align_rate.tsv) |
        mlr --itsv --omd cat >> visuliaztion/${PREFIX}_align_info.md
    echo -e "\n" >> visuliaztion/${PREFIX}_align_info.md

    cat remove_align/${PREFIX}/${DIR}/common.tsv >> visuliaztion/${PREFIX}_align_info.md
    echo -e "\n" >> visuliaztion/${PREFIX}_align_info.md

    head remove_bacteria/${DIR}/output/${PREFIX}/align_statistics.tsv |
      mlr --itsv --omd cat >> visuliaztion/${PREFIX}_align_info.md
    echo -e "\n" >> visuliaztion/${PREFIX}_align_info.md

    echo '![]'"(./plot/${PREFIX}/${PREFIX}_${DIR}.png)" >> visuliaztion/${PREFIX}_align_info.md
    echo -e "\n" >> visuliaztion/${PREFIX}_align_info.md
    echo '![]'"(./plot/${PREFIX}/${PREFIX}_${DIR}_CG.png)" >> visuliaztion/${PREFIX}_align_info.md
    echo -e "\n" >> visuliaztion/${PREFIX}_align_info.md
done

# rRNA_tRNA mRNA
for DIR in rRNA_tRNA mRNA;do
    echo -e "## ${DIR}\n" >> visuliaztion/${PREFIX}_align_info.md
    (head -n 1 visuliaztion/${PREFIX}_align_rate.tsv && grep "^${DIR}" visuliaztion/${PREFIX}_align_rate.tsv) |
        mlr --itsv --omd cat >> visuliaztion/${PREFIX}_align_info.md
    echo -e "\n" >> visuliaztion/${PREFIX}_align_info.md

    cat remove_align/${PREFIX}/${DIR}/common.tsv >> visuliaztion/${PREFIX}_align_info.md
    echo -e "\n" >> visuliaztion/${PREFIX}_align_info.md

    head remove_bacteria/${DIR}/output/${PREFIX}/align_statistics.tsv |
      tsv-select -e 2 |
      mlr --itsv --omd cat >> visuliaztion/${PREFIX}_align_info.md
    echo -e "\n" >> visuliaztion/${PREFIX}_align_info.md

    echo '![]'"(./plot/${PREFIX}/${PREFIX}_${DIR}.png)" >> visuliaztion/${PREFIX}_align_info.md
    echo -e "\n" >> visuliaztion/${PREFIX}_align_info.md
    echo '![]'"(./plot/${PREFIX}/${PREFIX}_${DIR}_CG.png)" >> visuliaztion/${PREFIX}_align_info.md
    echo -e "\n" >> visuliaztion/${PREFIX}_align_info.md
done

pandoc -i visuliaztion/${PREFIX}_align_info.md -o visuliaztion/${PREFIX}_align_info.html
rm visuliaztion/${PREFIX}_align_info.md
rm visuliaztion/${PREFIX}_align_rate.tsv
