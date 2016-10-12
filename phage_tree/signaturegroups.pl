#!/usr/bin/perl -w

#    Copyright 2001, 20002 Rob Edwards
#    For updates, more information, or to discuss the scripts
#    please contact Rob Edwards at redwards@utmem.edu or via http://www.salmonella.org/
#
#    This file is part of The Phage Proteome Scripts developed by Rob Edwards.
#
#    Tnese scripts are free software; you can redistribute and/or modify
#    them under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    They are distributed in the hope that they will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    in the file (COPYING) along with these scripts; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#


use DBI;
use strict;
open OUT, ">matches.html" || die "can't open matches.html for writing\n";

my $dbh=DBI->connect("DBI:mysql:phage", "apache") or die "Can't connect to database\n";

print OUT "<html><head><title>Phage group matches</title></head><body bgcolor=\"#FFFFFF\">\n";

my $tag; my %matches;
while (<>) {
  chomp;
  if (/^group/i) {
    if (%matches) {printdata(\%matches)}
     $tag = $_;
     undef %matches;
   }
 else {
   my @line = split /\s+/;
   my $file = shift @line;
   @line = sort {$a <=> $b} @line;
   my $newline = join "\t", @line;
   $matches{$newline}=1;
   }
  }
printdata(\%matches);
print OUT "</table>\n";
print OUT "</body></html>\n";



&niceexit(0);


sub niceexit {
	my $reason = shift;
	$dbh->disconnect;
	if ($reason) {print STDERR $reason; exit(-1)}
	else {exit(0)}
	}

sub printdata {
	my $matches = shift;
	my %matches = %$matches;
	print OUT "</table>\n<table border=1 width=\"100%\"><tr><td colspan=4><center><b>$tag</b></center></td></tr>\n";
      foreach my $key (keys %matches) {
   	my @line=split /\t/, $key;
	print OUT "<tr><td colspan=6><center><b>Signature proteins</b></center></td></tr>\n";
	print OUT "<tr><td>Bacteria</td><td>NCBI protein</td><td>Function</td><td>Protein sequence</td><td>Alignment</td><td>Protdist</td></tr>\n";

	foreach my $part (@line) {
   
      		my $exc = $dbh->prepare("select organism, proteinid, function, count from protein where count = '$part'") or &niceexit("$dbh->errstr");
		$exc->execute or &niceexit("$dbh->errstr");
		my @retrieved = $exc->fetchrow_array;
		unless ($retrieved[0]) {&niceexit("Can't find $part in protein\n")}

		$exc = $dbh->prepare("select organism from phage where count = '$retrieved[0]'") or &niceexit("$dbh->errstr");
		$exc->execute or &niceexit("$dbh->errstr");
		my @ret = $exc->fetchrow_array;
		unless ($ret[0]) {&niceexit("Can't find $retrieved[0] in phage\n")}
		my $orgcount=$retrieved[0];
		$retrieved[0] = $ret[0];
		my $count = pop(@retrieved);
		if ($retrieved[1]) {
		$retrieved[1]="\n" . '<a href="http://www.ncbi.nlm.nih.gov/entrez/viewer.cgi?val='.$retrieved[1].'">'.$retrieved[1].'</a>';
		}
		else {$retrieved[1]="&nbsp;"}
		print OUT "<tr><td>", join ("</td><td>", @retrieved);
		print OUT "<td><a href=\"/cgi-bin/phage.cgi?dbs=protein&count=$count&accession=1&locus=1&length=1&protein=1&protgenome=1&protseq=1\">$count</a>&nbsp;";
		print OUT '</td><td><a href="alignments/', $count,"_",$orgcount,'">',$count,"_",$orgcount,"</a></td>\n";
                print OUT '<td><a href="protdist/', $count,"_",$orgcount,'">',$count,"_",$orgcount,"</a></td>\n";
		print OUT "</tr>\n";
	}
	
	}
      print OUT "\n";
     }
    
