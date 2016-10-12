#!/usr/bin/perl -w

use strict;
use Getopt::Std;
my %opts;
getopts('o:f:q:', \%opts);

unless ($opts{f} && $opts{q} && $opts{o}) {
	die <<EOF;
$0
-f input file (fasta)
-q output quality file
-o output fastq file

EOF
}

print STDERR "Reading ... \n";
my $qu = read_fasta($opts{q}, 1);
my $fa = read_fasta($opts{f});
print STDERR "Done\n";

open(FQ, ">$opts{o}") || die "Can't open $opts{o} for writing";
foreach my $id (keys %$fa) {
	my $qual;
	foreach my $n (split /\s+/, $qu->{$id}) {
		$n += 33;
		my $c = chr($n);
		$qual.=$c;
	}
	print FQ "\@$id\n", $fa->{$id}, "\n+", "$id\n$qual\n";
}
close FQ;
		

=head1 read_fasta

Read a fasta format file and return a hash with the data, but this is not just limited to fasta sequence files - it will also handle quality scores and such.

Takes a file as argument, and an optional boolean to ensure that it is a quality file

usage: 
	my $fa=$rob->read_fasta($fastafile);
	my $fa=$rob->read_fasta("fastafile.gz"); # can also handle gzipped files
	my $qu=$rob->read_fasta($qualfile, 1);

Returns a reference to a hash - keys are IDs values are values

=cut

sub read_fasta {
	my ($file, $qual)=@_;
	if ($file =~ /.gz$/) {open(IN, "gunzip -c $file|") || die "Can't open a pipe to $file"}
	else {open (IN, $file) || die "Can't open $file"}
	my $f; my $t; my $s; my $newlinewarning;
	while (<IN>) {
		if (/\r/) 
		{
			print STDERR "The fasta file $file contains \\r new lines. It should not do this. Please complain bitterly.\n" unless ($newlinewarning);
			$newlinewarning=1;
			s/\r/\n/g; 
			s/\n\n/\n/g;
		}
		chomp;
		if (/^>/) {
			s#>##;
			if ($t) {
				if ($qual) {$s =~ s/\s+/ /g; $s =~ s/\s$/ /; $s =~ s/^\s+//}
				else {$s =~ s/\s+//g}
				$f->{$t}=$s;
				undef $s;
			}
			$t=$_;
		}
		else {$s .= " " . $_ . " "}
	}
	if ($qual) {$s =~ s/\s+/ /g; $s =~ s/\s$/ /; $s =~ s/^\s+//}
	else {$s =~ s/\s+//g}
	$f->{$t}=$s;
	close IN;
	return $f;
}


