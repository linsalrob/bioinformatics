#perl
#
# create a new genome_id.map file
#

use strict;

open(IN, "protein_locus_name.tsv") || die $!;
my %d; # d for data
while (<IN>) {
	chomp;
	my @a=split /\t/;
	$d{$a[0]}=\@a;
}
close IN;

my %g;
open(IN, "proteins_id.map") || die $!;
while (<IN>) {
	chomp;
	s/gbk\|//;
	s/refseq\|//;
	my @a=split /\t/;
	unless (defined $d{$a[1]}) {print STDERR "No data for $a[1]\n"; next}
	my ($prot, $genome) = split /\_/, $a[0];
	my $tg = $d{$a[1]}->[1];
	unless (defined $g{$genome}) {
		$g{$genome}=$tg;
		print join("\t", $genome, $tg), "\n";
	}
	if ($g{$genome} && $g{$genome} ne $tg) {
		print STDERR "Had $g{$genome} for $genome, but now have $tg\n";
	}
}

