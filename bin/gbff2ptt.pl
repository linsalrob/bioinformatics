#!/usr/bin/perl -w

# get protein information out of a genbank file

use Bio::SeqIO;
use strict;

my $usage=<<EOF;
$0 <list of genbankfiles>

EOF

die $usage unless ($ARGV[0]);

open(FA, ">fasta") || die "Can't open fasta";
open(PR, ">function") || die "Can't open function";

my $c;
foreach my $file (@ARGV)
{
	my $sio=Bio::SeqIO->new(-file=>$file, -format=>'genbank');
	while (my $seq=$sio->next_seq) {
		my $seqname=$seq->display_name;
		print STDERR "Parsing $seqname\n";
		foreach my $feature ($seq->top_SeqFeatures()) {
			my $id; # what we will call the sequence
			my ($trans, $gi, $geneid, $prod, $locus, $np, @xids);

			eval {$trans = join " ", $feature->each_tag_value("translation")};
			eval {$np = join " ", $feature->each_tag_value("protein_id")};
			next unless ($np);
			$c++;

			eval {
				foreach my $xr ($feature->each_tag_value("db_xref")) 
				{
					($xr =~ /GI/) ? ($gi = $xr) : 
					($xr =~ /GeneID/) ? ($geneid = $xr) : 
					(push @xids, $xr);
				}
			};

			eval {$locus = join " ", $feature->each_tag_value("locus_tag")};
			eval {$prod  = join " ", $feature->each_tag_value("product")};

			unless ($prod)  {print STDERR "No product for $np\n"; $prod="hypothetical protein"}
			my $strand = "";
			($feature->strand() > 0) ? ($strand = "+") : 
				($feature->strand() < 0) ? ($strand = "-") : 1;
			
			print FA ">$c\n$trans\n";
			print PR join("\t", $c, "refseq|$np", $feature->start(), $feature->end(), $strand,  $feature->location->to_FTstring(), 
			length($trans), $locus, $gi, $geneid, join(";", @xids), $prod), "\n";
		}
	}
}


