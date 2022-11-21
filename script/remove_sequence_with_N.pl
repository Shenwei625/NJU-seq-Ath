#! usr/bin/env perl
use strict;
use warnings;
use autodie;

=head1 NAME
remove_sequence_with_N.pl   --Remove sequence with N
=head1 SYNOPSIS
    perl remove_sequence_with_N.pl file.fa
=cut

open( my $in_fh1, "<", $ARGV[0] );

my @ALL;
while (<$in_fh1>) {
    chomp;
    push @ALL, $_;
}
close $in_fh1;

my $i = 0;
my @discard_sequence;
while ( $i <= $#ALL ) {
    if ( $ALL[$i] =~ /^>/ ) {
        my $sequence_name = $ALL[$i];
        $i++;
        my @all_seq;
        until ( ( $i > $#ALL ) or ( $ALL[$i] =~ /^>/ ) ) {
            my @seq = split( //, $ALL[$i] );
            foreach (@seq) {
                push @all_seq, $_;
            }
            $i++;
        }
        if ( "@all_seq" =~ m/N/ ) {
#            print "DISCARD $sequence_name\n";
            push @discard_sequence, $sequence_name;
        }
    }
}

my $j = 0;
while ( $j <= $#ALL ) { 
    if ( $ALL[$j] =~ /^>/ ) {
        my $sequence_name = $ALL[$j];
        $j++;
        my @all_seq2;
        until ( ( $j > $#ALL ) or ( $ALL[$j] =~ /^>/ ) ) {
            push @all_seq2, $ALL[$j];
            $j++;
        }
        if ( ! ( grep { $sequence_name eq $_ } @discard_sequence )) {
            print "$sequence_name\n";
            foreach (@all_seq2) {
                print "$_\n";
            }
        }
    }
}



