#!/usr/bin/perl -w

use strict;
use DBI;
use Bio::SearchIO;

my $dbh = DBI->connect('DBI:mysql:phage:n0000', 'rob', 'forestry') or die "Can't connect to database\n";

my ($indir, $outdir)=@ARGV;
unless ($indir && $outdir) {die "$0 <input directory of blasts> <output directory>"}

die "$outdir is not a directory" unless (-d $outdir);
die "$indir is not a directory" unless (-d $indir);
$indir =~ s/\/$//;
$outdir =~ s/\/$//;


my $outfile=$indir; $outfile =~ s/^.*\///;

my %details;
open(OUT, ">$outdir/$outfile") || die "can't opne $outdir/$outfile";
opendir(DIR, $indir) || die "Can't open $indir";
foreach my $file (grep {$_ !~ /^\./} readdir(DIR))
{
 my $sio=Bio::SearchIO->new(-file=>"$indir/$file", -format=>"blast");
 while (my $res=$sio->next_result)
 { 
  my $id = $res->query_name();
  while (my $hit=$res->next_hit)
  {
   my $hitname=$hit->name;
   my ($protein, $org);
   if ($hitname =~ /\_/) {($protein, $org)=split /\_/, $hitname}
   elsif ($hitname =~ /^\d+$/) {$org=$hitname}
   else {print STDERR "NOT SURE WHAT ORGANISM TO PARSE FROM $hitname\n"; next}
   unless ($details{$org})
   {
    my $exc=$dbh->prepare("select accession, locus, family, organism, count, genbankname from phage where count=$org") || die $dbh->errstr;
    $exc->execute|| die $dbh->errstr;
    my @ret=$exc->fetchrow_array;
    splice(@ret, 2, 0, "", "phage");
    $details{$org}=\@ret;
   }
   print OUT join("\t", $id, @{$details{$org}}, $hit->raw_score, $hit->significance), "\n";
  }
 }
} 

