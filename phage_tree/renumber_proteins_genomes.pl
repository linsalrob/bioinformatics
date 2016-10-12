#!/usr/bin/perl -w
#

use strict;
use Rob;
my $rob=new Rob;
use Getopt::Std;
my %opts = (
	'g' => 0,
	'p' => 0,
	'n' => 'proteins.function',
	'f' => 'proteins.fasta',
	's' => 'seqences.fasta',
	'o' => 'proteins_renumbered.fasta',
	'i' => 'proteins_id.map',
	'l' => 'phage_protein_genome_length.txt',
);

getopts('g:p:n:f:o:i:hs:l:', \%opts);

if ($opts{h}) {
	print <<EOF;
$0

Renumber a fasta file into the format we need for the phage tree.

-g start the genome count at this number (so you can concatenate with other files) (default $opts{g})
-p start the protein count at this number (so you can concatenate with other files) (default $opts{p})
-n name of the proteins.functions file that has tuples of [protein id, function, locus tag, genome] (default $opts{n})
-s name of the DNA sequences file (default $opts{s})
-f protein fasta file (default $opts{f})
-o output file for the fasta (default $opts{o})
-i output file for the idmap (default $opts{i})
-l output file for the proteins_genomes_lengths file (default $opts{l})
-h print this help
EOF
exit();
}

my $length = read_seq_len($opts{s});




open(IN, $opts{n}) || die $!;
my $genome=$opts{g};
my $last;
my %count;
my %gen;
my %len;
while (<IN>) {
	chomp;
	my @a=split /\t/;
	unless (defined $last) {
		$last=$a[3];
	}
	if ($a[3] ne $last) {
		$last=$a[3];
		$genome++;
	}
	$count{$a[0]}=$genome;
	$gen{$genome}=$a[3];
	$len{$genome}=$length->{$a[3]}
	# print STDERR "Adding |$a[0]|\n";
}
close IN;

open(IN, $opts{f}) || die $!;
open(OUT, ">$opts{o}") || die $!;
open(ID, ">$opts{i}") || die $!;
open(LEN, ">$opts{l}") || die $!;
my $c=$opts{p};
while (<IN>) {
	if (s/^>//) {
		chomp;
		if (defined $count{$_}) {
			print ID "${c}_$count{$_}\t", $gen{$count{$_}}, "\n";
			print LEN join("\t", "${c}_$count{$_}", $count{$_}, $len{$_}), "\n";
			print OUT ">", $c++, "_", $count{$_}, "\n";
		}
		else {
			print STDERR "|$_| NOT FOUND!!\n";
			print ID "${c}_00\t$_\n";
			print OUT ">", $c++, "_00\n";
		}
	}
	else {
		print OUT;
	}
}



sub read_seq_len {
	my $f=shift;
	my %l;
	my $fa = $rob->read_fasta($f);
	map {$l{$_}=length($fa->{$_})} keys %$fa;
	return \%l;
}
