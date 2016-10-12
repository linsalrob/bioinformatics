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
$| =1;

use constant IMAGE_WIDTH => 900;
use constant LEFT_MARGIN => 200;

my $dbh=DBI->connect("DBI:mysql:phage", "SQLUSER", "SQLPASSWORD") or die "Can't connect to database\n";
my $file= shift || &niceexit("drawtras.pl <file of prot transfers> <max horizontal size (kb)>\n");
my $maxsize = shift || &niceexit("drawtras.pl <file of prot transfers> <max horizontal size (kb)>\n");
my $longestgenome = $maxsize * 1000;

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
my $contig_bar_width = 30;
my $vspace = 10;
my $l_margin = LEFT_MARGIN;
my $t_margin = 10;
my $bgcolor = 'white';

my $image_height = $vspace + (@transfer+3) * ($contig_bar_width + $vspace) + 25;
my $image = Image->new(IMAGE_WIDTH, $image_height)->fill('white');
my %templateFrame;

#$longestgenome=1;
#foreach my $alllength (keys %$genomelength) {
#	if ($$genomelength{$alllength} > $longestgenome) {$longestgenome=$$genomelength{$alllength}}
#	}



my $offset = ${l_margin}-5;

# draw the size marker to scale
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
print STDERR "Adding $genome: ", int($$genomelength{$genome}/1000), " kb\n";
	my $origin;
	
	if ($genome > 1) {
		my $oldgenome = $genome-1;
		$origin = ($templateFrame{$oldgenome}->coords)[3];
		$origin += $vspace;
	}
	else {$origin = $t_margin+50}
	#generate the frame for each genome and add the organism name
	# if the length of the genome is longer than $longestgenome we need to add more height.
	my $heightfactor =1;
	if ($$genomelength{$genome} > $longestgenome) {
		#note unless the genome length is an exact factor, we need to add one to the heightfactor
		unless ($$genomelength{$genome}%$longestgenome) {$heightfactor = $$genomelength{$genome}/$longestgenome}
		else {$heightfactor = (int($$genomelength{$genome}/$longestgenome))+1}
		}
	#now convert the height into pixels for the image
	$heightfactor = $heightfactor * $contig_bar_width;
	$templateFrame{$genome} = $image->Frame->place(-coords =>  ["${l_margin}p", "${origin}p", 0.95, ".+${heightfactor}p"]);
	unless ($$organismname{$genome}) {print STDERR "WARNING: No genome name for $genome\n"}
	
	# generate the name and length. If the name is too long we want to divide it onto two lines.
	my $name = int($$genomelength{$genome}/1000);
	$name = $$organismname{$genome}." (".$name."kb)";

	if (length($name) > 27) {
		$name =~ /(.{0,27})\s+(.*?)$/;
		my $name1 = $1; my $name2 = $2;
		$templateFrame{$genome}->string(-coords => ["-${offset}p", '0p'], -text => $name1, -color => 'navy');
		$templateFrame{$genome}->string(-coords => ["-${offset}p", '.5'], -text => $name2, -color => 'navy');
		}
	else {
	$templateFrame{$genome}->string(-coords => ["-${offset}p", '0p'], -text => $name, -color => 'navy');
	}
	
	
	#add the arrows for the ORFS
	foreach my $orf (@{$$orfs{$genome}}) {
		my $y = 0.1;
		my ($start, $stop) = ($$beg{$orf}, $$end{$orf});
		if ($start > $stop) {($stop, $start) = ($start, $stop)}
		while (($start > $longestgenome) || ($stop>$longestgenome)) {
			$start -= $longestgenome;
			$stop -= $longestgenome;
			$y += (1/(int($$genomelength{$genome}/$longestgenome)+1));
			}

		$start = $start/$longestgenome;
		$stop = $stop/$longestgenome;
		if ($$complement{$orf}) {
		$templateFrame{$genome}->line(-coords=>[$start, $y, $stop, $y], -fill=>'red', -arrow=>'first', -width=>'0.05');
		}
		else {
		$templateFrame{$genome}->line(-coords=>[$start, $y, $stop, $y], -fill=>'red', -arrow=>'last', -width=>'0.05');
		}
	}
	
	# %lineat is used to store all the y positions so we can draw lines there. The value of $lineat will
	# end up being the highest stop position.
	# $axis_posn is used to store the position where to draw the y axis at, because in some cases this could be < 0
	my %lineat;
	my $axis_posn=0;
	# add the chart for the transfer frequency
	foreach my $gene (1 .. $#{$transfer[$genome]}) {
		if ($transfer[$genome][$gene]) {
			# $y is the maximum position on the chart. This will be the bottom of the box
			my $y=0.9;
			if ($$genomelength{$genome} > $longestgenome) {
				$y=(1/(int($$genomelength{$genome}/$longestgenome)+1));
				}
			# calculate the position of the rectangle
			my ($start, $stop) = ($$beg{$gene}, $$end{$gene});
			# adjust if the start or stop are longer than the scale
			while (($start > $longestgenome) || ($stop>$longestgenome)) {
				$start -= $longestgenome;
				$stop -= $longestgenome;
				$y+=(1/(int($$genomelength{$genome}/$longestgenome)+1));
				}
			$start = $start/$longestgenome;
			$stop = $stop/$longestgenome;
			
			# save axis position for later.
			if ($start < $axis_posn) {$axis_posn=$start}
			if ($stop < $axis_posn) {$axis_posn=$stop}
			# calculate the percent transfer. The one that is transferred the most will be 100%
			# note if there will be more than one line, we have to scale it down.
			my $height;
			if ($$genomelength{$genome} > $longestgenome) {
			  $height = $y-(($transfer[$genome][$gene]/$max{$genome})*
			     ((1/(int($$genomelength{$genome}/$longestgenome)+1))*0.6))
			  }
			else {$height = $y-(($transfer[$genome][$gene]/$max{$genome})*0.6)}
			
			# draw the damn box
			$templateFrame{$genome}->filledRectangle($start, $height, $stop, $y, 'blue');
			unless (exists $lineat{$y}) {$lineat{$y} = $stop}
			if ($lineat{$y} < $stop) {$lineat{$y} = $stop}

			}
		}
	# add the box around the chart. Do this last so it is on top!
	
	unless (%lineat) {$lineat{0.9} = $$genomelength{$genome}/$longestgenome}
	
	if ($$genomelength{$genome} > $longestgenome) {
		foreach my $y (keys %lineat) {
		  $templateFrame{$genome}->line(-coords=>[$axis_posn, $y, $lineat{$y}, $y], -fill=>'black', -width=>'0.01');
		}
		$templateFrame{$genome}->line(-coords=>[$axis_posn, '0.1', $axis_posn, '1'], -fill=>'black', -width=>'0.01');
		$templateFrame{$genome}->line(-coords=>["-${offset}p", '1.05', '1.5', '1.05'], -fill=>'grey', -width=>'0.01');
	}
	else {
		foreach my $y (keys %lineat) {
		  $templateFrame{$genome}->line(-coords=>[$axis_posn, $y, $lineat{$y}, $y], -fill=>'black', -width=>'0.01');
		}
		$templateFrame{$genome}->line(-coords=>[$axis_posn, '0.2', $axis_posn, '0.9'], -fill=>'black', -width=>'0.01');
		$templateFrame{$genome}->line(-coords=>["-${offset}p", '1.2', '1.5', '1.2'], -fill=>'grey', -width=>'0.01');
	}
}


  open IMG, ">$file.scaled.$maxsize.png";
  print IMG $image->png;
  close IMG;

&niceexit("file written to $file.scaled.$maxsize.png");






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



