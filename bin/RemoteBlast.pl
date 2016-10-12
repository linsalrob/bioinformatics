#!/usr/bin/perl -w 
#
# run remote blast against some sequences. Hopefully this will blast against the entire nr

use strict;
use lib '/home/rob/perl';
use Bio::Tools::Run::RemoteBlast;

my $usage =<<EOF;

$0 [-n,p,x,tn,s,tx] [-0] file [blast options]
        -n      blastn          nt query, nt database
        -p      blastp          protein query, protein database
        -x      blastx          nt query, protein database
        -tn     tblastn         protein query, nt database
        -s      psitblastn      protein query, nt database psi-blast
        -tx     tblastx         nt query, nt database
        
        -0      skip responses that don't match anything

EOF

my $blast=shift || die $usage;
my $blastprogram;
if ($blast eq "-n") {$blastprogram="blastn"}
elsif ($blast eq "-p") {$blastprogram="blastp"}
elsif ($blast eq "-x") {$blastprogram="blastx"}
elsif ($blast eq "-s") {$blastprogram="psiblast"}
elsif ($blast eq "-tn") {$blastprogram="tblastn"}
elsif ($blast eq "-tx") {$blastprogram="tblastx"}
die $usage unless ($blastprogram); 

my $test=shift || die $usage;

my @files;
if (-T $test) {push @files, $test}
elsif (-d $test) {
 opendir(DIR, $test) || die "Can't open dir $test\n";
 while (my $file2=readdir(DIR)) {
  push (@files, "$test/$file2");
 }
}
else {
 die "Not sure what $test is";
} 

foreach my $file (@files) {
my $rb=Bio::Tools::Run::RemoteBlast->new(
  '-prog' => $blastprogram, '-data' => 'nr', '-expect' => 10, '-readmethod' => 'SearchIO'
);

 my $res = $rb->submit_blast($file); 
 
 my $matched=0;
 #sleep 5;
 #my @rids=$rb->each_rid;
 #until ($matched == scalar @rids) {
 while (my @rids=$rb->each_rid) {
  #@rids=$rb->each_rid;
  foreach my $rid (@rids) {
    print STDERR "For $file, there are ", scalar @rids, " rids waiting. $matched are complete\n";
    my $rc = $rb->retrieve_blast($rid);
   if (!ref($rc)) {
      if ($rc < 0) {
         $rb->remove_rid($rid);
      }
      sleep 5;
   } else {
    $matched++;
    $rb->save_output("$file.$blastprogram.nr");
    $rb->remove_rid($rid);
   }
  }
 }
}
