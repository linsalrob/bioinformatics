#!/usr/bin/perl -w

use strict;
use Getopt::Std;
my %opts;
getopts('i:f:q:', \%opts);

unless ($opts{i}) {
	die <<EOF;
$0
-i input file (fastq)
-f output fasta file (use - for stdout)
-q output quality file (use - for stderr) (optional)

The -f and -q are optional

Note: If you redirect the output to STDOUT you explicitly have to say if you want the quality put out to stderr. This is so you can pipe the output

EOF
}

my $fastqf;
my ($o, $q);
if ($opts{i} =~ /fastq$/) {
	$o = $q = $opts{i};
	$o =~ s/fastq/fasta/;
	$q =~ s/fastq/qual/;
	open($fastqf, $opts{i}) || die "can't open $opts{i}";
}
elsif ($opts{i} =~ /fastq.gz$/) {
	$o = $q = $opts{i};
	$o =~ s/fastq.gz/fasta/;
	$q =~ s/fastq.gz/qual/;
	open($fastqf, "gunzip -c $opts{i}|") || die "can't open pipe to gunzip $opts{i}";
}
elsif ($opts{i} =~ /fastq.zip$/) {
	$o = $q = $opts{i};
	$o =~ s/fastq.zip/fasta/;
	$q =~ s/fastq.zip/qual/;
	open($fastqf, "unzip -p $opts{i}|") || die "can't open pipe to unzip $opts{i}";
}
else {
	die "can't figure out $opts{i}";
}

$opts{f} ? 1 : ($opts{f} = $o);
if ($opts{f} ne "-") {$opts{q} ? 1 : ($opts{q} = $q)}


if ($opts{f} eq "-") {
	open(FA, ">&STDOUT") || die "Can't open a pipe to stdout";
} else {
	open(FA, ">$opts{f}") || die "Can't open $opts{f} for writing";
}

if (!$opts{q}) {
	open(QU, ">/dev/null") || die "can't open a pipe to /dev/null";
}
elsif ($opts{q} eq "-") {
	open(QU, ">&STDERR") || die "Can't open a pipe to stderr";
} else {
	open(QU, ">$opts{q}") || die "Can't open $opts{q} for writing";
}
my $l=0; my $n=-1;
while (<$fastqf>) {
	$n++;
	if ($n == 0) {
		s/^@/>/;
		print FA;
		print QU;
	}
	elsif ($n == 1) {
		print FA;
	}
	elsif ($n == 2) {
	}
	elsif ($n == 3) {
		chomp;
		foreach my $c (split //, $_) {
			my $x = ord($c);
			$x -= 33;
			print QU $x, " ";
		}
		print QU "\n";
		$n = -1
	}
}
close FA; close QU; close $fastqf;
