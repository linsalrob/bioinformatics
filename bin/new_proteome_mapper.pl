#!/usr/bin/perl -w

use strict;

unless ($ARGV[0]) {
	die "$0 [-m minimum peptide length]  <list of files>";
}

my @files;
my $min = 0;

while (@ARGV) {
	my $t=shift @ARGV;
	if ($t eq "-m") {$min=shift @ARGV}
	elsif (-e $t) {push @files, $t}
	else {print STDERR "rur ro! What is option $t?\n"}
}


# a new version of the proteome mapper. This will take the files from Eidy and compare them to the nr
# ARGV is a list of text files, and we'll map them to those ids

my $pept;
foreach my $file (@files) {
	my $id = $file;
	$id =~ s/\.txt//;
	open(IN, $file) || die "cant open $file";
	while (<IN>) {
		m/^(\S+)/;
		if ($1) {
			my $peptide = $1;
			$pept->{$peptide}->{$id}=1 if (length($peptide) >= $min);
		}
	}
	close IN;
}

# now process the nr
my ($seq, $id)=(undef, undef);
open(NR, "fastacmd -d nr -D 1 |") || die "can't open a pipe to fastacmd";
while (<NR>) {
	chomp;
	if (index($_, ">") == 0) {
		$seq && (&process($id, $seq)); 
		$seq = undef;
		$id = $_;
		$id =~ s/^\>//;
	}
	else {$seq .= $_}
}
$seq && (&process($id, $seq));


sub process {
	my ($id, $seq)=@_;
	foreach my $p (keys %$pept) {
		if (index($seq, $p) > -1) {print join("\t", $p, $id, join("; ", keys %{$pept->{$p}}), index($seq, $p), $seq), "\n"}
	}
}



