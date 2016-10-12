#!/usr/bin/perl -w

# read in output from Beltrans program and the ss_key file to convert back to subsystems

use strict;

my $usage=<<EOF;
$0
-k subsystems key file (default: ~/seed_database/ss_key.txt) (optional)
-g subsystems group file (default is not defined) (optional)
-i input file from Beltran


EOF

my ($keyf, $inf,$groupf)=('','','');
while (@ARGV) {
 my $t=shift;
 if ($t eq "-k") {$keyf=shift}
 elsif ($t eq "-i") {$inf=shift}
 elsif ($t eq "-g") {$groupf=shift}
}

$keyf = "/home/rob/seed_database/ss_key.txt" unless ($keyf);
die $usage unless ($inf && $keyf);

my $outf=$inf; 

($outf =~ s/\.txt$/\.named/) ? (1) : ($outf .= ".named");

my $group;
if (open(IN, $groupf)) {
 while (<IN>)  {
  chomp;
  my @a=split /\t/;
  $a[0] =~ s/\_/ /g; $a[0]=~s/  / /g; $a[0]=~s/^\s+//; $a[0]=~s/\s+$//;
  unless (defined $a[2]) {$a[2]=''}
  unless (defined $a[3]) {$a[3]=''}
  $group->{$a[0]}=$a[2]."\t".$a[3] unless (exists $group->{$a[0]});
 }
}

my $ss; 
open(IN, $keyf) || die "Can't open $keyf";
while (<IN>) {
 chomp;
 my @a=split /\t/;
 $ss->{$a[1]}=$a[0];
}
close IN;
my $conf;
my %output; my %allsamples;
open(TXT, ">$outf.txt") || die "Can't open $outf.txt for writing";
open (IN, $inf) || die "Can't open $inf";
while (<IN>) {
 if (/confidence/) {$conf=$_; print TXT; next}
 if (/^\s+(\d+)\s+(.*)$/) {
  my ($id, $sample)=($1, $2);
  $allsamples{$sample}++;
  if ($ss->{$id}) {
   my $output=$ss->{$id} . "\t" . $sample;

   # add subsystem group information
   my $a=$ss->{$id};  $a =~ s/\_/ /g; $a=~s/  / /g; $a=~s/^\s+//; $a=~s/\s+$//;
   if ($group->{$a}) {
    next if ($group->{$a} =~ "Delete"); 
    $output = $group->{$a} ."\t". $output;
   }
   else {
    if ($groupf) {
     print STDERR "No group found for $a looking in |$groupf|\n";
     $output = "Unclassified\t\t".$output;
    }
   }
   push (@{$output{$sample}}, $output);
  }
  else {print TXT "$id\t$sample\n"; print STDERR "No subsystem from $_\n"}
 }
 else {print TXT}
}

print TXT map {join"\n", @{$output{$_}}, ''} sort keys %output;


open(HTM, ">$outf.html") || die "Can't open $outf.html for writing";
my $samplehtml=join " and ", map {"$_ ($allsamples{$_} subsystems)"} keys %allsamples;
print HTM <<EOF;
<html>
<head>
<title>Comparative Analysis</title>
</head>
<body>

<p>Comparison of $samplehtml with $conf</p>

<table border=1>

EOF

if ($groupf) {print HTM "<tr><th colspan=2>Classification</th><th>Subsystem</th><th>Sample</th></tr>\n"} 
else {print HTM "<tr><th>Subsystem</th><th>Sample</th></tr>\n"}


my @colors=('#98FB98', '#FFB6C1', '#DEB887');
my %color;
map {$color{$_}=shift @colors} keys %allsamples;

foreach my $sam (keys %output) {
 foreach my $line (sort @{$output{$sam}}) {
  $line =~ s/\_/ /g; $line =~ s/  / /g;
  $line =~ /.*\t(.*?)$/;
  my $color=$color{$1};
  $line =~ s#\t#</td><td>#g;
  print HTM "<tr style='background-color: $color'><td>$line</td></tr>\n";
 }
}

print HTM "</table>\n</html>\n";

