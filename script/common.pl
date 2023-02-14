#! usr/bin/env perl
use strict;
use warnings;
use autodie;
use Getopt::Long;

=head1 NAME
    common.pl   --Output the same lines in two files
=head1 SYNOPSIS
    perl common.pl --index in.tsv --search in2.tsv
        Options:
            --help\-h   Brief help message
            --index\-i  The smaller of the two files
            --search\-s The larger of the two files
            --stdin     Get content from STDIN. It will not been not valid with a provided '--search'
=cut

Getopt::Long::GetOptions(
    'help|h'     => sub { Getopt::Long::HelpMessage(0) },
    'index|i=s' => \my $index,
    'search|s=s'   => \my $search,
    'stdin'      => \my $stdin,
) or Getopt::Long::HelpMessage(1);

open( my $in_fh1, "<", $index );
my %ALL;
while (<$in_fh1>) {
    $ALL{$_} = 1;
}
close $in_fh1;

if ( defined($search) ) {
    open( my $in_fh2, "<", $search );
    while (<$in_fh2>) {
        print "$_" if exists $ALL{$_};
    }
    close $in_fh2;
}
elsif ( defined($stdin) ) {
    while (<>) {
        print "$_" if exists $ALL{$_};
    }
}
else {
    die("You should provide stdin!");
}

__END__

