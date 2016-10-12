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
use OGD;

use constant IMAGE_WIDTH => 900;
use constant LEFT_MARGIN => 200;


my $dbh=DBI->connect("DBI:mysql:phage", "SQLUSER", "SQLPASSWORD") or die "Can't connect to database\n";
my $file= shift || &niceexit("drawtras.pl <file of prot transfers>\n");

my @transfer; my %max;
open (IN, "$file") || &niceexit("Can't open $file\n");
while (<IN>) {
	chomp;
	my ($genome, $gene, $count) =split /\t/;
	$transfer[$genome][$gene] = $count;
	unless ($max{$genome}) {$max{$genome}=$transfer[$genome][$gene]}
	if ($transfer[$genome][$gene] > $max{$genome}) {$max{$genome}=$transfer[$genome][$gene]}
}
close IN;

my ($genomelength, $beg, $end, $organismname, $orfs, $complement) = &getdata(\%max);

#generate the image. This is largely taken from projection_map.cgi (part of ERGO).
my $template_color = 'tan';
my $template_highlight_color = 'wheat';
my $contig_bar_width = 25;
my $vspace = 10;
my $l_margin = LEFT_MARGIN;
my $t_margin = 10;
my $bgcolor = 'white';

my $image_height = $vspace + (@transfer+3) * ($contig_bar_width + $vspace) + 25;
my $image = Image->new(IMAGE_WIDTH, $image_height)->fill('white');
my %templateFrame;

my $longestgenome = 1;
foreach my $alllength (keys %$genomelength) {
	if ($$genomelength{$alllength} > $longestgenome) {$longestgenome=$$genomelength{$alllength}}
	}

my $offset = ${l_margin}-5;

{
$templateFrame{'sizemarker'} = $image->Frame->place(-coords =>  ["${l_margin}p", "${t_margin}p", 0.95, ".+${contig_bar_width}p"]); 
my $fivekb = 5000/$longestgenome;
$templateFrame{'sizemarker'}->line(-coords=>['0', '1', $fivekb, '1'], -fill=>'black', -width=>'0.01');
$templateFrame{'sizemarker'}->line(-coords=>['0', '0.75', '0', '1'], -fill=>'black', -width=>'0.01');
$templateFrame{'sizemarker'}->line(-coords=>[$fivekb, '0.75', $fivekb, '1'], -fill=>'black', -width=>'0.01');
$templateFrame{'sizemarker'}->string(-coords => ['0', '0'], -text => '0', -color => 'navy');
$templateFrame{'sizemarker'}->string(-coords => [$fivekb, '0'], -text => '5 kb', -color => 'navy');
$templateFrame{'sizemarker'}->line(-coords=>["-${offset}p", '1.2', '1.5', '1.2'], -fill=>'grey', -width=>'0.01');
}

foreach my $genome (1 .. $#transfer) {
print STDERR "Adding $genome\n";
	my $origin;
	
	if ($genome > 1) {
		my $oldgenome = $genome-1;
		$origin = ($templateFrame{$oldgenome}->coords)[3];
		$origin += $vspace;
	}
	else {$origin = $t_margin+50}
	#generate the frame for each genome and add the organism name
	$templateFrame{$genome} = $image->Frame->place(-coords =>  ["${l_margin}p", "${origin}p", 0.95, ".+${contig_bar_width}p"]);
	unless ($$organismname{$genome}) {print STDERR "WARNING: No genome name for $genome\n"}
	
	my $name = int($$genomelength{$genome}/1000);
	$name = $$organismname{$genome}." (".$name."kb)";

	if (length($name) > 27) {
		$name =~ /(.{0,27})\s+(.*?)$/;
		my $name1 = $1; my $name2 = $2;
		print STDERR "Name1 : $name1 Name2: $name2\n";
		$templateFrame{$genome}->string(-coords => ["-${offset}p", '0p'], -text => $name1, -color => 'navy');
		$templateFrame{$genome}->string(-coords => ["-${offset}p", '.5'], -text => $name2, -color => 'navy');
		}
	else {
	$templateFrame{$genome}->string(-coords => ["-${offset}p", '0p'], -text => $name, -color => 'navy');
	}
	# add a line around the genome, and a line for the transfer charts;
	
	#add the arrows for the ORFS
	foreach my $orf (@{$$orfs{$genome}}) {
		my ($start, $stop) = ($$beg{$orf}, $$end{$orf});
		if ($start > $stop) {($stop, $start) = ($start, $stop)}
#		$start = $start/$$genomelength{$genome};
#		$stop = $stop/$$genomelength{$genome};
		$start = $start/$longestgenome;
		$stop = $stop/$longestgenome;
		if ($$complement{$orf}) {
		$templateFrame{$genome}->line(-coords=>[$start, '0', $stop, '0'], -fill=>'red', -arrow=>'first', -width=>'0.05');
		}
		else {
		$templateFrame{$genome}->line(-coords=>[$start, '0', $stop, '0'], -fill=>'red', -arrow=>'last', -width=>'0.05');
		}
	}	
	# add the proteins and their arrows
	foreach my $gene (1 .. $#{$transfer[$genome]}) {
		if ($transfer[$genome][$gene]) {
			# calculate the percent transfer. The one that is transferred the most will be 100%
			my $height = 1-(($transfer[$genome][$gene]/$max{$genome})*0.8);
			# now add the rectangles for the transfer frequency
			my ($start, $stop) = ($$beg{$gene}, $$end{$gene});
#			$start = $start/$$genomelength{$genome};
#			$stop = $stop/$$genomelength{$genome};
			$start = $start/$longestgenome;
			$stop = $stop/$longestgenome;
		
			$templateFrame{$genome}->filledRectangle($start, $height, $stop, '1', 'blue');

			}
		}
	my $genomefraction = $$genomelength{$genome}/$longestgenome;
	$templateFrame{$genome}->line(-coords=>['0', '1', "$genomefraction", '1'], -fill=>'black', -width=>'0.01');
	$templateFrame{$genome}->line(-coords=>['0', '0.2', '0', '1'], -fill=>'black', -width=>'0.01');
	$templateFrame{$genome}->line(-coords=>["-${offset}p", '1.2', '1.5', '1.2'], -fill=>'grey', -width=>'0.01');
	
	}


  open IMG, ">$file.scaled.png";
  print IMG $image->png;
  close IMG;

&niceexit();






sub niceexit {
	my $reason = shift;
	$dbh->disconnect;
	if ($reason) {print STDERR $reason; exit(-1)}
	else {exit(0)}
	}

sub getdata {
print STDERR "Getting data\n";
	my $genomes = shift;
	my %genomelength; my %start; my %stop; my %organismname; my %orfs; my %complement;
	my $exc = $dbh->prepare("SELECT count, organism, sequence from phage" ) or croak $dbh->errstr;
	$exc->execute or die $dbh->errstr;
	while (my @retrieved = $exc->fetchrow_array) {
			$organismname{$retrieved[0]} = $retrieved[1];
			$genomelength{$retrieved[0]} = length($retrieved[2]);
			}
	
	
	foreach my $key (keys %$genomes) {
print STDERR "\tGetting $key\n";
		$exc = $dbh->prepare("SELECT count,start,stop,complement from protein where organism = $key" ) or croak $dbh->errstr;
		$exc->execute or die $dbh->errstr;
		while (my @retrieved = $exc->fetchrow_array) {
			($start{$retrieved[0]}, $stop{$retrieved[0]}) = ($retrieved[1], $retrieved[2]);
			push (@{$orfs{$key}}, $retrieved[0]);
			$complement{$retrieved[0]} = $retrieved[3];
			}
		}
print STDERR "\tDone\n";
	return \%genomelength, \%start, \%stop, \%organismname, \%orfs, \%complement;
}



