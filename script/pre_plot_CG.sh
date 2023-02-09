#! usr/bin/bash
# bash pre_plot.sh PREFIX
PREFIX=$1

echo "===> RAW"
J=../data/${PREFIX}
pigz -dcf $J/R1.fq.gz | perl -e'
  while (<>) {
    chomp( my $seq_name = $_ );
    chomp( my $seq = <> );
    my $seq_length = length($seq);
    my $C_count = $seq =~ tr/C/C/;
    my $G_count = $seq =~ tr/G/G/;
    my $CG_count = ( ( ( $C_count + $G_count ) / $seq_length ) * 100 );
    chomp( my $info = <> );
    chomp( my $quality = <> );
    printf "%.0f\n", $CG_count;
  }
' | sort -n | uniq -c | perl -ne'
      s/^\s+//;
      print "$_";
    ' | tr " " "\t" | tsv-select --fields 2,1 > $J/CG_distribution.tsv

(echo -e "CG_content(%)\tTotal_number" && cat $J/CG_distribution.tsv) > tem&&
  mv tem $J/CG_distribution.tsv

for DIR in rRNA_conserve rRNA_tRNA mRNA;do
echo "===> ${DIR}"
pigz -dcf ${DIR}/filter/${PREFIX}/R1* | perl -e'
  while (<>) {
    chomp( my $seq_name = $_ );
    chomp( my $seq = <> );
    my $seq_length = length($seq);
    my $C_count = $seq =~ tr/C/C/;
    my $G_count = $seq =~ tr/G/G/;
    my $CG_count = ( ( ( $C_count + $G_count ) / $seq_length ) * 100 );
    chomp( my $info = <> );
    chomp( my $quality = <> );
    printf "%.0f\n", $CG_count;
  }
' | sort -n | uniq -c | perl -ne'
      s/^\s+//;
      print "$_";
    ' | tr " " "\t" | tsv-select --fields 2,1 > ${DIR}/filter/${PREFIX}/CG_distribution.tsv

(echo -e "CG_content(%)\tKeep_number" && cat ${DIR}/filter/${PREFIX}/CG_distribution.tsv) > tem&&
  mv tem ${DIR}/filter/${PREFIX}/CG_distribution.tsv
done

echo "===> merge"
perl script/tsv_join_plus.pl ../data/${PREFIX}/CG_distribution.tsv rRNA_conserve/filter/${PREFIX}/CG_distribution.tsv > rRNA_conserve/filter/${PREFIX}/keep_remove_CG.tsv
(echo -e "CG_content(%)\tTotal_number\tKeep_number\tKeep_proportion" && sed '1d' rRNA_conserve/filter/${PREFIX}/keep_remove_CG.tsv) > tem&&
    mv tem rRNA_conserve/filter/${PREFIX}/keep_remove_CG.tsv

perl script/tsv_join_plus.pl rRNA_conserve/filter/${PREFIX}/CG_distribution.tsv rRNA_tRNA/filter/${PREFIX}/CG_distribution.tsv > rRNA_tRNA/filter/${PREFIX}/keep_remove_CG.tsv
(echo -e "CG_content(%)\tTotal_number\tKeep_number\tKeep_proportion" && sed '1d' rRNA_tRNA/filter/${PREFIX}/keep_remove_CG.tsv) > tem&&
    mv tem rRNA_tRNA/filter/${PREFIX}/keep_remove_CG.tsv

perl script/tsv_join_plus.pl rRNA_tRNA/filter/${PREFIX}/CG_distribution.tsv mRNA/filter/${PREFIX}/CG_distribution.tsv > mRNA/filter/${PREFIX}/keep_remove_CG.tsv
(echo -e "CG_content(%)\tTotal_number\tKeep_number\tKeep_proportion" && sed '1d' mRNA/filter/${PREFIX}/keep_remove_CG.tsv) > tem&&
    mv tem mRNA/filter/${PREFIX}/keep_remove_CG.tsv

echo "===> pre_plot"
for DIR in rRNA_conserve rRNA_tRNA mRNA;do
echo -e "CG_content\tnumber\tgroup" > ${DIR}/filter/${PREFIX}/plot_CG.tsv

sed '1d' ${DIR}/filter/${PREFIX}/keep_remove_CG.tsv | perl -ne'
  chomp;
  if (/^(\S+)\t(\S+)\t(\S+)/) {
    my $reads_length = $1;
    my $total_number = $2;
    my $keep_number = $3;
    my $remove_number = ( $total_number - $keep_number );
    print "$reads_length\t$keep_number\tKeep\n";
    print "$reads_length\t$remove_number\tRemove\n";
  }
' >> ${DIR}/filter/${PREFIX}/plot_CG.tsv

tsv-filter -H --ge CG_content:25 --le CG_content:75 ${DIR}/filter/${PREFIX}/plot_CG.tsv > tem&&
    mv tem ${DIR}/filter/${PREFIX}/plot_CG.tsv
done
