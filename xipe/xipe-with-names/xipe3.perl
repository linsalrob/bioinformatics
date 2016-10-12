#!/usr/bin/perl -s

#Copyright (C) 2005 Beltran Rodriguez-Brito,
#Pat MacNarnie, and Rober Edwards
#
#This program is free software; you can redistribute it and/or
#modify it under the terms of the GNU General Public License
#as published by the Free Software Foundation; either version 2
#of the License, or (at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program; if not, write to the Free Software
#Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#
#Released as part of the xipe package as supplemental online material
#to "A statistical approach for metagenome comparisons"
# by Beltran Rodriguez-Brito, Forest Rohwer, and Robert A. Edwards.

use strict;

if(!@ARGV){
	print "Usage: ./slow.perl mediansfile rangesfile repeats confidence name1 name2\n";
	exit(); }
my $filename  = $ARGV[0]; my $filename2 = $ARGV[1]; 
my $therepeats = $ARGV[2]; my $confidence = $ARGV[3];
my $name1 = $ARGV[4]; my $name2 = $ARGV[5];
my $j; my $k; my $thelength; my $keepgoing; my $confi; my $xleft; my $xright;

open(FD, "<$filename") or die("Couldn't open mediansfile\n") ;
open(FH, "<$filename2") or die("Couldn't open rangesfile\n") ;
my @arr; my $maria; my $petra; my $yleft; my $mymedian; my $yright;
$j = 1;
$keepgoing = 1;
print $confi;
print " $confidence confidence\n";
while( ($maria = <FD>)&&($petra=<FH>)&&($keepgoing==1) ) {
        $maria =~ s/\n|\r//g;
        $petra =~ s/\n|\r//g;
	my @lines = split(/ /, $maria);
        my @linesrange = split(/ /, $petra);
        $thelength = @lines;
        my $vergatmp = int( $thelength/2);
        if ( $thelength == $therepeats ) {
        $xright = int( $confidence*$thelength + 0.5 );
        $xleft = int($thelength-$xright );
        $yleft = @linesrange[$xleft];
        $yright = @linesrange[$xright];
        $mymedian = @lines[ $vergatmp ];
        #print " left,$yleft,median,$mymedian,right,$yright \n ";
        if ( $mymedian < $yleft ) {
             print " $j $name1\n"; 
             }
        if ( $mymedian > $yright ) {
             print " $j $name2\n";
             }

#        if ( @lines[ int($thelength/2) ] < @linesrange[$xleft] ) {
#             print " 1  $j ,  \n"; 
#             }
#        if ( @lines[ int($thelength/2) ] > @linesrange[$xright] ) {
#             print " 2 ", $j , " \n";
#             }
        $j = $j + 1;
        }
}
print " Done!\n";
close(FD);
close(FH);
