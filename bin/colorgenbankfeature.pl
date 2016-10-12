#!/usr/bin/perl -w

# color a single feature in a genbank file

use strict;
use Bio::SeqIO;
use Getopt::Std;

my %opts;
getopts('s:d:', \%opts);
my $usage = <<EOF;
$0 
-s source directory
-d destination directory
EOF


die $usage unless ($opts{s} && $opts{d});


foreach my $file (`find $opts{s} -name \*gbk`) {
	$file =~ /^(.*)\/(.*?gbk)$/;
	my ($dir, $ofile)=($1, $2);
	unless ($dir && $ofile) {die "couldn't parse directory/file from $file"}
	
	$dir=~s/$opts{s}/$opts{d}/;
	unless (-e "$dir") {`mkdir -p $dir`}

	my $sio=Bio::SeqIO->new(-file=>$file, -format=>'genbank');
	my $sout = Bio::SeqIO->new(-file=>">$dir/$ofile", -format=>'genbank');
	while (my $seq=$sio->next_seq) {
		my $seqname=$seq->display_name;

		print STDERR "Parsing $seqname\n";
		foreach my $feature ($seq->top_SeqFeatures()) {
			# check and see if has color
			my $color;
			eval {$color = join " ", $feature->each_tag_value("color")};
			next if ($color && $color !~ /^\s+$/);

			my $prod;
			eval {$prod  = join " ", $feature->each_tag_value("product")};
			next unless (defined $prod);

			if (($prod =~ /site-specific/i && $prod =~ /recombinase/i) || ($prod =~ /integrase/i)) {
				$feature->add_tag_value('color', '2');
				print STDERR "Added 2 to $file\n";
			}
			elsif ($prod =~ /phage/i || $prod =~ /tail/i || $prod =~ /capsid/i || $prod =~ /portal protein/i) {
				$feature->add_tag_value('color', '3');
				print STDERR "Added 3 to $file\n";
			}
		}
		$sout->write_seq($seq);
	}
}


