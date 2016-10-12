#!/usr/bin/perl -w

use strict;
use DBI;

my $dbh=DBI->connect("DBI:mysql:phage", "rob", "forestry") or die "Can't connect to database\n";

foreach my $type (qw[Streptococcus Staphylococcus vibrio campylobacter listeria]) {
 print "*** $type ***\n";
 my $exc=$dbh->prepare("select count, organism, host from phage where (organism like '\%$type\%' or host like '\%$type\%')") || die $dbh->errstr;
 $exc->execute || die $dbh->errstr;
 while (my @res=$exc->fetchrow_array) {
  print join "\t", @res, "\n";
 }
}

