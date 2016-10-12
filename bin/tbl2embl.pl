#!/usr/bin/perl -w
#

use strict;
use Bio::SeqIO;
use Bio::SeqFeature::Generic;
use Bio::Seq;
use Getopt::Std;

my %opts;
getopts('f:t:o:', \%opts);

unless ($opts{f} && $opts{t} && $opts{o}) {
	die <<EOF;
$0
-f fasta file (required)
-t five column tab-separated genbank file (required). This is the output from gbf2tbl.pl
-o output file (required). This will be in EMBL format.

EOF
}

# read the fasta file and store that as a hash
my $seqin = Bio::SeqIO->new(-file=>$opts{f}, -format=>"fasta");
my %fasta;
while (my $s = $seqin->next_seq()) {
	$fasta{$s->display_id} = $s->seq;
}

# create the output file
my $seqout = Bio::SeqIO->new(-file=>">$opts{o}", -format=>"embl");

# sequences that we have included in the output
my %seen;

# read the table and store all of this as a hash
my $id = undef;
my $currfeat = undef;
my $seqobj;
open(IN, $opts{t}) || die "can't open $opts{t}";
while (<IN>) {
	chomp;
	if (/^>Feature\s+(\S+.*)$/) {
		my $newid = $1;
		if ($id) {
			$currfeat && $seqobj->add_SeqFeature($currfeat);
			$seqout->write_seq($seqobj);
		}
		$id = $newid;
		if (!$fasta{$id}) {
			die "We don't have a sequence for $id in the fasta file at $_\n";
		}

		$currfeat = undef;
		$seqobj = Bio::Seq->new(-display_id => $id, -seq => $fasta{$id});
		$seen{$id}=1;

		next;
	}

	unless ($id) {die "Don't have an ID to continue with. Does your tab separated table contain a line that begins '>Feature '?\n"}

	my @parts = split("\t");
	if ($parts[0] && $parts[1] && $parts[2]) {
		if ($currfeat) {
			$seqobj->add_SeqFeature($currfeat);
		}
		my $strand = 1;
		if ($parts[1] < $parts[0]) {($parts[0], $parts[1], $strand) = ($parts[1], $parts[0], -1)}
		$currfeat = new Bio::SeqFeature::Generic(
			-start       => $parts[0],
			-end         => $parts[1],
			-strand      => $strand,
			-primary_tag => $parts[2],
		);
	}
	else {
		$currfeat->add_tag_value($parts[3], $parts[4]);
	}
}
$currfeat && $seqobj->add_SeqFeature($currfeat);
$seqout->write_seq($seqobj);



foreach  my $id (keys %fasta) {
	unless ($seen{$id}) {
		print STDERR "We didn't find any features for $id\n";
		$seqobj = Bio::Seq->new(-display_id => $id, -seq => $fasta{$id});
		$seqout->write_seq($seqobj);
	}
}


