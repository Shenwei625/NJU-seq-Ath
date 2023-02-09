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
    chomp( my $info = <> );
    chomp( my $quality = <> );
    print "$seq_length\n";
  }
' | sort -n | uniq -c | perl -ne'
      s/^\s+//;
      print "$_";
    ' | tr " " "\t" | tsv-select --fields 2,1 > $J/length_distribution.tsv
  
(echo -e "reads_length\tTotal_number" && cat $J/length_distribution.tsv) > tem&&
  mv tem $J/length_distribution.tsv

for DIR in rRNA_conserve rRNA_tRNA mRNA;do
echo "===> ${DIR}"
pigz -dcf ${DIR}/filter/${PREFIX}/R1* | perl -e'
  while (<>) {
    chomp( my $seq_name = $_ );
    chomp( my $seq = <> );
    my $seq_length = length($seq);
    chomp( my $info = <> );
    chomp( my $quality = <> );
    print "$seq_length\n";
  }
' | sort -n | uniq -c | perl -ne'
      s/^\s+//;
      print "$_";
    ' | tr " " "\t" | tsv-select --fields 2,1 > ${DIR}/filter/${PREFIX}/length_distribution.tsv

(echo -e "reads_length\tKeep_number" && cat ${DIR}/filter/${PREFIX}/length_distribution.tsv) > tem&&
  mv tem ${DIR}/filter/${PREFIX}/length_distribution.tsv
done

echo "===> merge"
perl script/tsv_join_plus.pl ../data/${PREFIX}/length_distribution.tsv rRNA_conserve/filter/${PREFIX}/length_distribution.tsv > rRNA_conserve/filter/${PREFIX}/keep_remove.tsv
perl script/tsv_join_plus.pl rRNA_conserve/filter/${PREFIX}/length_distribution.tsv rRNA_tRNA/filter/${PREFIX}/length_distribution.tsv > rRNA_tRNA/filter/${PREFIX}/keep_remove.tsv
perl script/tsv_join_plus.pl rRNA_tRNA/filter/${PREFIX}/length_distribution.tsv mRNA/filter/${PREFIX}/length_distribution.tsv > mRNA/filter/${PREFIX}/keep_remove.tsv

echo "===> pre_plot"
for DIR in rRNA_conserve rRNA_tRNA mRNA;do
echo -e "reads_length\tnumber\tgroup" > ${DIR}/filter/${PREFIX}/plot.tsv

sed '1d' ${DIR}/filter/${PREFIX}/keep_remove.tsv | perl -ne'
  chomp;
  if (/^(\S+)\t(\S+)\t(\S+)/) {
    my $reads_length = $1;
    my $total_number = $2;
    my $keep_number = $3;
    my $remove_number = ( $total_number - $keep_number );
    print "$reads_length\t$keep_number\tKeep\n";
    print "$reads_length\t$remove_number\tRemove\n";
  }
' >> ${DIR}/filter/${PREFIX}/plot.tsv

head -n 83 ${DIR}/filter/${PREFIX}/plot.tsv > tem&&
  mv tem ${DIR}/filter/${PREFIX}/plot.tsv
done
