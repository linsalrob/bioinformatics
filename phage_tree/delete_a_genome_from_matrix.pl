#!/usr/bin/perl -w
#

use strict;

=pod

delete a single genome from a matrix. This is to be used, for example, with matrix or matrix.nusubreplicates and assumes the matrix is in phylip format

The required arguments are the name of the matrix file and the number of the genome to be deleted. This should be the row number that the genome is in, where the first row in the file is row 0

the optional argument is the output file. If no file is provided, stdout will be used.

=cut

use Getopt::Std;
my %opts = ('o' => '-');
getopts('m:n:o:', \%opts);

unless ($opts{m} && $opts{n}) {
	die <<EOF;
	$0
-m matrix file
-n line number of genomes to ignore. This is 0 based with the first line of the file (the line with the number of genomes on it) being zero
-o output file (optional : if not provided STDOUT will be used

EOF

}


my $first = 1;
my $row=0;
open(IN, $opts{m}) || die "Can't open $opts{m}";
open(OUT , ">$opts{o}") || die "can't open $opts{o}";
while (<IN>) {
	chomp;
	if ($first) {
		$_--;
		print OUT "$_\n";
		$first=undef;
		next;
	}
	$row++;
	s/^(\S+\s+)//;
	my $id=$1;
	my @a=split /\s+/;
	if ($row == $opts{n}) {
		print STDERR "Deleting $id. Before : ", $a[$opts{n}-2], " Splicing : ", $a[$opts{n}-1], " After: ", $a[$opts{n}], "\n";
	}
	my @omit = splice @a, $opts{n}-1, 1;
	if ($row == $opts{n}) {
		print STDERR "Deleting $id. Omit: @omit\n";
		next;
	}
	
	print OUT $id, join(" ", @a), "\n";
}
	






