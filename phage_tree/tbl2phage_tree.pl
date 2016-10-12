#!/usr/bin/perl -w
#

use strict;
use Rob;
my $rob = new Rob;
use Getopt::Std;
my %opts = (
	'd' => 'sequences.fasta',
	'p' => 'proteins.fasta',
	'n' => 'proteins.function',
	'l' => 'phage_protein_genome_length.txt',
);

getopts('d:p:n:l:h', \%opts);

if ($opts{h}) {die <<EOF;
$0

-d DNA sequence file (default $opts{d})
-p protein fasta file (default $opts{p})
-n protein functions file (default $opts{n})
-h print this help menu
EOF
}

my $len = read_len($opts{d});

# create phage_protein_genome_length

open(IN, $opts{n}) || die "Can'topen $opts{n}";
open(OUT, ">$opts{l}") || die "Can't write to $opts{l}";
while (<IN>) {
	chomp;
	my @a=split /\t/;
	print OUT join("\t", $a[0], $a[3], $len->{$a[3]}), "\n";
}
close IN;
close OUT;





sub read_len {
	my $f = shift;
	my $fa = $rob->read_fasta($f);
	my %l;
	map {$l{$_}=length($fa->{$_})} keys %$fa;
	return \%l;
}


