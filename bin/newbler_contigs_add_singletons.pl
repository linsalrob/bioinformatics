#!/usr/bin/perl -w

# combine the all contigs and the assembly file. 

use strict;
use Getopt::Std;
use Rob;

my %opts;
getopts('f:a:o:', \%opts);

my $usage = <<EOF;
$0
-f fasta file of original sequences
-a all contigs file from 454
-o file to write the output to
-r 454ReadStatus.txt file

EOF

die $usage unless ($opts{f} && $opts{a} && $opts{o});
unless ($opts{r}) {
	$opts{r} = $opts{a};
	$opts{r} =~ s/454.*?Contigs.fna/454ReadStatus.txt/;
	if ($opts{a} eq $opts{r} || !-e $opts{r}) {
		die "You didn't specify a valid 454ReadStatus.txt file. We tried $opts{r}\n$usage\n";
	}
}
my %single;

my %singletonReasons = (
	"Singleton" => 1,
	"TooShort"  => 1,
	"Outlier"   => 1,
);

open(IN, $opts{r}) || die "can't open $opts{r}";
while (<IN>) {
	chomp;
	my @a=split /\t/;
	if ($singletonReasons{$a[1]}) {
		$single{$a[0]}=1;
	}
	elsif ($#a < 2) {
		print STDERR "Had too few columns, so kept $a[0] as a singleton and added '$a[1]' as a reason to keep singletons\n";
		$singletonReasons{$a[1]}=1;
		$single{$a[0]}=1;
	}
}

close IN;

my $fa = Rob->read_fasta($opts{f});
map {my $x=$_; s/\s+.*$//; $fa->{$_}=$fa->{$x}} keys %$fa;

open(OUT, ">$opts{o}") || die "Can't write to $opts{o}";
foreach my $id (keys %single) {
	unless ($fa->{$id}) {print STDERR "WARNING NO SEQUENCE FOR $id\n"; next}
	print OUT ">$id\n", $fa->{$id}, "\n";
}
open(IN, "$opts{a}") || die "Can't open $opts{a}";
print OUT while (<IN>);
close IN;
close OUT;

