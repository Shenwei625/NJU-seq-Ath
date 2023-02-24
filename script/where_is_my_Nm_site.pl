#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

=head1 NAME
    where_is_my_Nm_site.pl   --output the position of the Nm site on the gene
=head1 SYNOPSIS
    perl where_is_my_Nm_site.pl ath.gff3 Ath_${TISSUE}_mrna_Nm_score.tsv
=cut


my @ALL;
open( my $info_fh, "<",  $ARGV[0]);
while (<$info_fh>) {
    chomp;
    push @ALL, $_;
}
close($info_fh);

my %gene_info;
my $i = 0;
while ( $i < $#ALL ) {
    if ($ALL[$i] =~ /ID=gene:(\w+);/) {
        my $gene_name = $1;
        $i++;
        until ( ($ALL[$i] =~ /ID=gene:(\w+);/) or ( $i >= $#ALL )) {
            my @info = split( /\t/, $ALL[$i] );
            my $transcript_id = $1 if $info[8] =~ /transcript:((\w+\.?[0-9]*))/;
            push ( @{ $gene_info{$gene_name} }, $info[2]);
            push ( @{ $gene_info{$gene_name} }, $info[3]);
            push ( @{ $gene_info{$gene_name} }, $info[4]);
            push ( @{ $gene_info{$gene_name} }, $transcript_id);
            $i++;
        }
        if ( $i == $#ALL ) {
            my @info2 = split( /\t/, $ALL[$i] );
            my $transcript_id = $1 if $info2[8] =~ /transcript:((\w+\.?[0-9]*))/;
            push ( @{ $gene_info{$gene_name} }, $info2[2]);
            push ( @{ $gene_info{$gene_name} }, $info2[3]);
            push ( @{ $gene_info{$gene_name} }, $info2[4]);
            push ( @{ $gene_info{$gene_name} }, $transcript_id);
        }
    }
}

# foreach my $trans ( sort( keys( %gene_info ) ) ) {
#    print "$trans\n";
# }

# find location
open( my $sc_fh, "<",  $ARGV[1]);
print "Gene\tsite\ttranscript\tlocation\tstart\tend\n";
while (<$sc_fh>) {
    my @info3 = split( /\t/, $_ );
    my $loc_gene = $info3[6];
    my $A=0;
    my $B=1;
    my $C=2;
    my $D=3;
    while ($C < $#{ $gene_info{$loc_gene} }) {
    if ( ( $info3[1] >= @{ $gene_info{$loc_gene} }[$B] )  and ( $info3[1] <= @{ $gene_info{$loc_gene} }[$C] ) ) {
        print "$info3[6]\t$info3[1]\t@{ $gene_info{$loc_gene} }[$D]\t@{ $gene_info{$loc_gene} }[$A]\t@{ $gene_info{$loc_gene} }[$B]\t@{ $gene_info{$loc_gene} }[$C]\n";
        $A = $A + 4;
        $B = $B + 4;
        $C = $C + 4;
        $D = $D + 4;
    } 
    else {
        $A = $A + 4;
        $B = $B + 4;
        $C = $C + 4;
        $D = $D + 4;
    }
    }
}

close($sc_fh);

__END__
