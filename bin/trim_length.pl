#!/usr/bin/perl -w

use strict;
use Getopt::Std;
use Rob;
my $rob=new Rob;

my %opts;
getopts('q:f:x:n:o:', \%opts);

die <<EOF unless (($opts{q} || $opts{f}) && ($opts{x} || $opts{n}));
$0
-f fasta file OR
-q fastq file
-x maximum length of sequence to include
-n minimum length of sequence to include
-o output file (or else STDOUT)

EOF

unless (defined $opts{n}) {$opts{n} = 0}
unless (defined $opts{x}) {$opts{x} = 1e20}
unless (defined $opts{o}) {$opts{o}="-"} # print to SDTOUT
open(OUT, ">$opts{o}") || die "can't open $opts{o}";

if ($opts{f}) {
	my $fa=$rob->read_fasta($opts{f});
	map {print OUT ">$_\n", $fa->{$_}, "\n" if (length($fa->{$_}) > $opts{n} && length($fa->{$_}) < $opts{x})} keys %$fa;
	close OUT;
}
if ($opts{q}) {
	my $tple = $rob->read_fastq($opts{q});
	map {print OUT "\@$_\n", $tple->{$_}->[0], "\n+$_\n", $tple->{$_}->[1], "\n" if (length($tple->{$_}->[0]) > $opts{n} && length($tple->{$_}->[0]) < $opts{x})} keys %$tple;
}
