#!/usr/bin/env perl 

# find repeats in all of the phage genomes.

# the repeats must be within 2kb of the start and end. 
# Actually: 
#   left_start must be 1,
#   left_end must be Sherwood phage start + 100
#   right_start must be sherwood phage end - 100
#   right_end must be length of sequence

use strict;
use Data::Dumper;
use lib '/home3/redwards/bioinformatics/Modules';
use RepeatFinder;
use Bio::SeqIO;

use Getopt::Std;
my %opts;
getopts('f:o:r:i:j:k:l:', \%opts);

unless ($opts{f} && $opts{o}) {
	die <<EOF;
$0
-f sequence file (filetype will be guessed)
-o output file
-r length of repeats (default=10)
-i left end start for searching (default = 0)
-j left end stop  for searching (default = length of the sequence)
-k right end start for searching (default = 0)
-l right end stop for searching (default = length of the sequence)
EOF
}



my $seqin=Bio::SeqIO->new(-file=>$opts{f});
my $seqout=Bio::SeqIO->new(-file=>">$opts{o}", -format=>'genbank');


my %args;
if ($opts{r}) {$args{"-minimum"}=$opts{r}} else {$args{"-minimum"}=10}
if ($opts{i}) {$args{"-left_start"}=$opts{i}}
if ($opts{j}) {$args{"-left_end"}=$opts{j}}
if ($opts{k}) {$args{"-right_start"}=$opts{k}}
if ($opts{l}) {$args{"-right_end"}=$opts{l}}

while (my $seq=$seqin->next_seq) {
	$args{"-seq"}=$seq;

	my $rep=RepeatFinder->new(%args);
	my $newseq=$rep->joined_repeats;
	$seqout->write_seq($newseq);
}



