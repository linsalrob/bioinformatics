#!/usr/bin/perl -w
#

use strict;
use lib '/home3/redwards/bioinformatics/Modules/';
use Rob;
my $rob=new Rob;
use Getopt::Std;
my %opts;
getopts('l:r:n:o:p:', \%opts);

my $usage = <<EOF;
-l left reads (file that ends in 1)
-r right reads (file that ends in 2)

-n number of reads to sample

-o left reads output file
-p right reads output file

All parameters are required.

EOF

die $usage unless ($opts{l} && $opts{r} && $opts{n} && $opts{o} && $opts{p});


print STDERR "Reading $opts{l}\n";
my $f = $rob->read_fastq($opts{l});
print STDERR "Reading $opts{r}\n";
my $g = $rob->read_fastq($opts{r});
print STDERR "Shuffling\n";
my $keys = $rob->rand([keys %$f]);
print STDERR "Writing\n";
open(L, ">$opts{o}") || die "$! $opts{o}";
open(R, ">$opts{p}") || die "$! $opts{p}";

foreach my $k (splice(@$keys, 0, $opts{n})) {
	my $t = $k;
	$t =~ s/1$/2/;
	print L "\@$k\n", $f->{$k}->[0], "\n+\n", $f->{$k}->[1], "\n";
	print R "\@$t\n", $g->{$t}->[0], "\n+\n", $g->{$t}->[1], "\n";
}
	


