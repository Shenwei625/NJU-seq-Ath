#! usr/bin/perl
use strict;
use warnings;
use autodie;

my @motif_len = qw(0 138 138 138 138 98 109 92 86 70 59);

open( my $fh1, "<", $ARGV[0] );

while (<$fh1>) {
    my @info = split( /\t/, $_ );
    my $site = 1;
    my @len_info = split( /-/, $info[1] );
    foreach (@len_info) {
        if (/^\[\+(\d+)\]/) {
            print "$info[0]:";
            print "$site-";
            $site = ( $site + $motif_len[$1] -1 );
            print "$site\n";
            $site++;
        }elsif (/^(\d+)/) {
            $site = ( $site + $1 );
        }
    }
}

__END__
