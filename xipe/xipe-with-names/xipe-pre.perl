#! /usr/bin/perl


# xipe-pre.pl
# Florent Angly
# Converts Human-readable data to Xipe-redable data


use warnings;
use strict;


# Usage
my $usage =<<STRING;
$0 converts human-readable data to Xipe-readable data (at the subsystem and role level).
Usage: $0 -i blasfile1 -j blastfile2 -k keyfile OPTIONS
Options:
   -f filterfile: keep only the PEGs specified in the file (one per line)
STRING


# Get arguments
my ($filetoconvert1, $filetoconvert2, $keyfile, $filterfile);
if (scalar @ARGV == 0) {die "$usage"};
while (@ARGV) {
   my $arg = shift @ARGV;
   if ("$arg" eq "-i") {
      $filetoconvert1= shift @ARGV;
   } elsif ("$arg" eq "-j") {
      $filetoconvert2 = shift @ARGV;
   } elsif ("$arg" eq "-k") {
      $keyfile = shift @ARGV;
   } elsif ("$arg" eq "-f") {
      $filterfile = shift @ARGV;
   } else {
      die "$usage";
   }
}
my $tempsub = "$keyfile.$filetoconvert1.$filetoconvert2.sub";
my $temprole = "$keyfile.$filetoconvert1.$filetoconvert2.role";


# Make sure the presence of the mandatory arguments are enforced
if (not($filetoconvert1) || not($filetoconvert1) || not($keyfile)) {
   die "$usage";
}


# Load fig vs sub vs role hashes
my (%subsystemhash, %rolehash);
open FIG, $keyfile || die "Cannot open file '$keyfile': $!\n\n";
while (<FIG>) {
   if (m/^(.*)\t(.*)\t(.*)$/) {
      #  SUBSYSTEM translation
      $subsystemhash{$1} = $3;
      # ROLE translation
      $rolehash{$1} = $2;
   } else {
      die "Cannot parse file '$keyfile': Bad file format\n\n";
   }
}
close FIG;


# Load PEGs from filter file if necessary
my @filterpegs;
if ($filterfile) {
   open(FILTERFILE, "$filterfile") || die("Cannot open file \"$filterfile\": $!");
   @filterpegs = <FILTERFILE>; 
}


# Main comp
my @tmpsub = ();
my @tmprole = ();
my %corresp; # Key is Fig ID, value 1 is subsystem name, value 2 is subsystem number, value 3 is role name, and value 4 is role number
my @filelist = ($filetoconvert1, $filetoconvert2);
my $subnumber = '1';
my $rolenumber = '2';
my $seqid = '';
my $oldseqid = '';
my $maxsubnumber = $subnumber;
my $maxrolenumber = $rolenumber;
foreach my $filein(@filelist) {
   my @results_sub = ();
   my @results_role = ();
   open FILEIN, $filein || die "Cannot open file '$filein': $!\n\n";
   while (<FILEIN>) {
      my $line = $_; 
      if ($line =~ m/^(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)$/) {

         $seqid = $1;
         my $fig = $2;
         my $subname = $subsystemhash{$2};
         my $rolename = $rolehash{$2};

         # Take only the best hits
         if (not ($seqid eq $oldseqid)) {
            
            # Filter only the lines mentioned in the filterfile if a filterfile was specified
            my $linetokeep;
            if ($filterfile) { # if there is a filterfile
               $linetokeep = 0; # do not keep any line by default
               foreach my $figtofilter(@filterpegs) { 
                  chomp $figtofilter;
                  if ( $figtofilter eq $fig) {
                     $linetokeep = 1; # keep this line
                  }
                  #print "$fig\t$figtofilter\t$linetokeep\n";
               }  
            } else { # if there is no filterfile
               $linetokeep = 1; #keep all lines
            }



            #
            if ($linetokeep == 1) { 
		            # Check the ouputfile for the existence of subsystem
		            #print "$seqid";
		            my $subfound = 0; # not found
		            foreach my $entry(@tmpsub) {
		               if ($entry =~ m/^(.*)\t(.*)$/) {
		                  my $oldsubname = $1;
		                  my $oldsubnumber = $2;
		                  # Record maximum subsystem number
		                  if ($oldsubnumber > $maxsubnumber) {
		                     $maxsubnumber = $oldsubnumber;
		                  }
		                  # if the subsystem name is known, take the same subsystem number
		                  if ($oldsubname eq $subname) {
		                     $subnumber = $oldsubnumber;
		                     $subfound = 1;
		                     last;
		                  }
		               }
		            }
		            # otherwise invent a new number for the subsystem and save it
		            if ($subfound == 0) {
		               $subnumber = $maxsubnumber + 1;
		               my $newline = "$subname\t$subnumber\n";
		               push @tmpsub, $newline; # Add to array 
		            }
		                        
		            # Check the ouputfile for the existence of role
		            #print "$seqid";
		            my $rolefound = 0; # not found
		            foreach my $entry(@tmprole) {
		               if ($entry =~ m/^(.*)\t(.*)$/) {
		                  my $oldrolename = $1;
		                  my $oldrolenumber = $2;
		                  # Record maximum role number
		                  if ($oldrolenumber > $maxrolenumber) {
		                     $maxrolenumber = $oldrolenumber;
		                  }
		                  # if the role name is known, take the same role number
		                  if ($oldrolename eq $rolename) {
		                     $rolenumber = $oldrolenumber;
		                     $rolefound = 1;
		                     last;
		                  }
		               }
		            }
		            # otherwise invent a new number for the role and save it
		            if ($rolefound == 0) {
		               $rolenumber = $maxrolenumber + 1;
		               my $newline = "$rolename\t$rolenumber\n";
		               push @tmprole, $newline; # Add to array 
		            }
		            # New result line (input for Xipe) is:
		            push @results_sub, "$seqid\t$subnumber\t$rolenumber\n"; # subsystems are in the second column
		            push @results_role, "$seqid\t$rolenumber\t$subnumber\n"; # roles are in the second column
		            #print "$seqid\t$subnumber\t$rolenumber\n";
            }
         }
         $oldseqid = $seqid;
      } else {
         die "Cannot parse input file '$filein': Bad file format\n\n";
      }
   }
   close FILEIN;
   # Save result arrays in a file
   open SUBOUT, ">$filein.sub.pre" || die "Cannot write file '$filein.sub.pre': $!\n\n";
   foreach my $line(@results_sub) {
      print SUBOUT "$line";
   }
   close SUBOUT;
   open ROLEOUT, ">$filein.role.pre" || die "Cannot write file '$filein.role.pre': $!\n\n";
   foreach my $line(@results_role) {
      print ROLEOUT "$line";
   }
   close ROLEOUT;
}



open(TMPSUB, ">$tempsub") || die "Cannot write file '$tempsub': $!\n\n";
foreach my $line(@tmpsub) {print TMPSUB "$line"};
close TMPSUB;

open(TMPROLE, ">$temprole") || die "Cannot write file '$temprole': $!\n\n";
foreach my $line(@tmprole) {print TMPROLE "$line"};
close TMPROLE;


exit;
