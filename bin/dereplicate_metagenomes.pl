#!/usr/bin/perl -w

use strict;
use Getopt::Std;

my %opts;
getopts('f:q:', \%opts);

unless ($opts{f}) {
	die "$0\n-f fasta file\n-q quality file\n";
}

my $filename;
if ($opts{f} =~ /^(.*)\.fna$/) {$filename=$1}
elsif ($opts{f} =~ /^(.*)\.fa$/) {$filename=$1}
elsif ($opts{f} =~ /^(.*)\.fasta$/) {$filename=$1}
else {die "Can't figure out whether $opts{f} is a fasta file (it doesn't end fna, fa, or fasta)"}


# read the fasta file to find out if sequneces are substrs
my ($duplicate, $five, $three)=(0,0,0);
open(THREE, ">three_prime_duplicates.fasta") || die "cant write to 3' dups";
open(FIVE, ">five_prime_duplicates.fasta") || die "cant write to 3' dups";
open(EXACT, ">exact_duplicates.fasta") || die "can;t write to exact dups";
my $seqs; my $allseqs;
my %fives; my %threes;
open(IN, $opts{f}) || die "can't open $opts{f}";
my $seq; my $tag;
my %ignore;
while (<IN>) {
	chomp;
	if (/^>/) {
		&clean($seq, $tag);
		$tag=$_;
		$seq="";
	} else {
		$seq .= uc($_);
	}

}
&clean($seq, $tag);
print STDERR ":After dereplicating there were $duplicate duplicate sequences, $five 5' overlapping sequences, and $three 3' overlapping sequences\n";

my $qu;
if ($opts{q}) {
	open(IN, $opts{q}) || die "can't open $opts{q}";
	while (<IN>) {
		chomp;
		if ($ignore{$_}) {$tag=""; next}
		if (/^\>/) {$tag=$_; next}
		if ($tag) {$qu->{$tag} .= $_}
	}
}


open(FA, ">$filename.dereplicated.fa") || die "can't open $filename.dereplicated.fa";
if ($qu) {open(QU, ">$filename.dereplicated.qu") || die "Can't open $filename.dereplicated.qu";}
foreach my $k (keys %$seqs) {
	next if ($ignore{$k});
	if ($qu && $qu->{$k}) {
		print QU "$k\n",$qu->{$k},"\n";
	}
	print FA "$k\n",$seqs->{$k},"\n";
}




sub clean {
	my ($seq, $tag)=@_;
	if ($seq) {
		if ($allseqs->{$seq}) {
		$duplicate++; 
		$ignore{$tag}=1;
		print EXACT $allseqs->{$seq}, "\n$seq\n$tag\n$seq\n\n";
}
		else {
			my $fp = substr($seq, 0, 20);
			my $tp = substr($seq, -20);
			if ($fives{$fp}) {
				foreach my $key (@{$fives{$fp}}) {
					if (!$seqs->{$key}) {die "No sequence $key"}
					if (index($seqs->{$key}, $seq) >= 0) {
# $seqs->{$key} is longer
						print FIVE "$key\n", $seqs->{$key}, "\n$tag\n$seq\n\n";
						$five++;
						$ignore{$tag}=1;
					}
					elsif (index($seq, $seqs->{$key}) >= 0) {
# $seq is longer
						print FIVE "$key\n", $seqs->{$key}, "\n$tag\n$seq\n\n";
						$five++;
						$ignore{$key}=1;
					}
				}
			}
			push @{$fives{$fp}}, $tag;

			if ($threes{$tp}) {
				foreach my $key (@{$threes{$tp}}) {
					if (index($seqs->{$key}, $seq) >= 0) {
# $seqs->{$key} is longer
						print THREE "$key\n", $seqs->{$key}, "\n$tag\n$seq\n\n";
						$three++;
						$ignore{$tag}=1;
					}
					elsif (index($seq, $seqs->{$key}) >= 0) {
# $seq is longer
						print THREE "$key\n", $seqs->{$key}, "\n$tag\n$seq\n\n";
						$three++;
						$ignore{$key}=1;
					}
				}
			}
			push @{$threes{$tp}}, $tag;
		}
			$seqs->{$tag}=$seq;
			$allseqs->{$seq}=$tag;
	}
}


