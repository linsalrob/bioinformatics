#!/usr/bin/perl -w

# run a deconvlution program to remove duplicates from all the sequences, and then run an assembly program to 
# assemble all the data.

# arguments: a fasta file, a quality scores file, and a directory to put the results in

use strict;

my $deconv = '/home/redwards/bioinformatics/bin/dereplicate_metagenomes.pl';
my $assembler = '/usr/local/454/version2.5/bin/runAssembly';

use Getopt::Std;
my %opts;
getopts('f:q:d:', \%opts);


unless ($opts{f}) {
	die <<EOF;
$0
-q quality file
-f fasta file
EOF
}
# first deconvolute the sequences, renaming the output

my $filename;
if ($opts{f} =~ /^(.*)\.fna$/) {$filename=$1}
elsif ($opts{f} =~ /^(.*)\.fa$/) {$filename=$1}
elsif ($opts{f} =~ /^(.*)\.fasta$/) {$filename=$1}
else {die "Can't figure out whether $opts{f} is a fasta file (it doesn't end fna, fa, or fasta)"}

if (!$opts{q}) {$opts{q}=""}
`$deconv -f $opts{f} -q $opts{q}`;
`$assembler $filename.dereplicated.fa`;




