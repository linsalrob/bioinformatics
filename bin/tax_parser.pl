#!/usr/bin/perl -w

# create a gi -> tax parser
# based on parse.cgi

# Also note the following rules to be nice to NCBI:
#       * Run retrieval scripts on weekends or between 9 PM and 5 AM ET weekdays for any series of more than 100 requests.
#       * Send E-utilities requests to http://eutils.ncbi.nlm.nih.gov, not the standard NCBI web address.
#       * Make no more than one request every 3 seconds.
#       * Use the URL parameter email, and tool for distributed software, so we can track your project and contact you if there is a problem.
#       * NCBI's Disclaimer and Copyright notice must be evident to users of your service.  NLM does not claim the copyright on the abstracts
#               in PubMed; however, journal publishers or authors may. NLM provides no legal advice concerning distribution of copyrighted 
#               materials, consult your legal counsel.



use strict;
use LWP::Simple;
use lib '/home/rob/perl';
use Bio::SearchIO;


my $usage="$0 <blast dir> OPTIONS\n\t-l use link\n\t-f use efetch\n\t-c convert GI to UI\n";

my $blastdir=shift || die $usage;
my $args = join " ", @ARGV;

my ($el, $ef)=(0,0);
if ($args =~ /-l/) {$el=1}
elsif ($args =~ /-f/) {$ef=1}

my $url;
my $qparam;
if ($el) 
{
 print STDERR "Using elink\n";
 $url='http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?';
 %$qparam=(
    'dbfrom'   => "nucleotide",
    'db'       => "taxonomy",
    'cmd'      => "neighbor",
    'field'    => '',
    'retmode'  => 'XML',
    'retstart' => '0',
    'retmax'   => '20',
    'tool'     => 'BLAST PARSER',
    'email'    => 'rob@salmonella.org'
 );

 #   'dbto'     => "taxonomy",   note: db works better!
 #   'cmd'      => "neighbor",


}
elsif ($ef)
{
 print STDERR "Using efetch\n";
 $url='http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?';
 %$qparam=(
   'db'       => "nucleotide",
   'field'    => '',
   'retmode'  => 'XML',
   'retstart' => '0',
   'rettype'  => 'gb',
   'tool'     => 'BLAST PARSER',
   'email'    => 'rob@salmonella.org'
 ); 
}
else
{
 die $usage;
}






my $time=0;
opendir(DIR, "$blastdir") || die "Can't open $blastdir";
while (my $file=readdir(DIR))
{
 next if ($file =~ /^\./);
 my $searchio=Bio::SearchIO->new(-file=>"$blastdir/$file", -format=>"blast");
 while (my $result = $searchio->next_result) 
 {
  my $query; my $posn=0;
  my %description; my %significance; my %raw_score;
  while (my $hit = $result->next_hit) 
  {
   #if ($hit->description =~ /$keyword/i) 
   # these are the things that we can store about the Hit:
   # name, description, length, algorithm, raw_score, significance, rank, 
   

   
   my @q=split /\|/, $hit->name;
   $$query[$posn]=$q[1];
   $description{$posn}=$hit->description;
   $raw_score{$posn}=$hit->raw_score;
   $significance{$posn}=$hit->significance;
   $posn++;
  }
  next unless ($query && scalar @$query);
  my $uid;
  if ($args =~ /-c/) {($time, $uid)=convert($query, $time)}
  #my $results=all_at_once($uid, $qparam);
  my $results=one_at_a_time($uid, $url, $qparam);
  my @results=split /\n/, $results;
  my ($id, $tax)=(0,0);
  my @id; my %taxid;
  foreach (@results)
  {
   if (m#<IdList>#) {$id=1; $tax=0; next}
   elsif (m#</IdList>#) {$id=0; $tax=0; next}
   elsif (m#<LinkSetDb>#) {$id=0; $tax=1; next}
   elsif (m#</LinkSetDb>#) {$id=0; $tax=0; next}
   next unless ($id || $tax);
   next unless (m#<Id>(\d+)</Id>#);
   my $found=$1;
   if ($id) {unshift @id, $found} 
   else 
   {
    my $match=shift @id;
    if ($taxid{$match}) {unless ($found eq $taxid{$match}) {print STDERR "Rewriting taxid for $match from $taxid{$match} to $found\n"}}
    $taxid{$match}=$found;
   }
  }

#print STDERR "taxids: \n", map {"$_=>$taxid{$_}\n"} keys %taxid;
  foreach my $i (0 .. $#$query)
  {
   print join "\t", $$query[$i], $$uid{$$query[$i]}, $taxid{$$uid{$$query[$i]}}, &tax_data($taxid{$$uid{$$query[$i]}}), $description{$i}, $raw_score{$i}, $significance{$i}, "\n";
  }
 }
}

exit(0);




sub convert {
 my ($gb, $time)=@_;
 # we need to convert gb to ui numbers
 # this is the command:
 # http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nucleotide&retmode=XML&rettype=uilist&term=AF121970.1
 my %ui;
 my @gb=@$gb; # note if we just shift from @$gb it really messes things up later!
 return 0 unless (scalar @$gb);
 my ($retstart, $retmax)=(0,$#$gb+5);
 my $convert="http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nucleotide&retmode=XML&retstart=$retstart&retmax=$retmax&rettype=uilist&term=".(join ",", @$gb);
 while (time-$time<3) {sleep 1} 
 my $result=LWP::Simple::get($convert);
 $time=time;
 #print "Trying to find $#$gb UI's from: $result\n";
 while ($result =~ s#^.*?<Id>(\d+)</Id>##s)
 {
  my $ui=$1;
  my $gi=shift @gb;
  $ui{$gi}=$ui;
 }
 #print "Found $#ui UI's\n";
 return ($time, \%ui);
}


sub all_at_once {
 # tHIS IS THE RIGHT WAY TO GET THE DATA, BUT UNFORTUNATELY ncbi
 # doesn't really support this, and the data sucks

  my ($uid, $url, $qparam)=@_; 
  $$qparam{'id'}=join ",", values %$uid;
  $$qparam{'retmax'}=(scalar keys %$uid)+5;
  my $link=$url .  join "&", map {"$_=$$qparam{$_}"} keys %$qparam;
  # wait until 3 seconds have passed since the last request. Shouldn't be too hard, but we should be nice
  # even though we are requesting 100 articles at a time
  while (time-$time<3) {sleep 1}
  $time=time;
  print STDERR "\n\nLINK\n\n$link\n\n";
  # there is a problem - the NCBI page doesn't return the taxid's in any meaningful form, and therefore we have to request them one at
  # a time

  
  my $results=LWP::Simple::get($link);
  delete $$qparam{'id'};
 return $results;
}


sub one_at_a_time {
 # this is a bad way to get the data. Oh well!
 # we're going to whack the crap out of NCBI for a short bit
 # they should return the data in a meaningful way.

 my $results; # all results
 # we'll do the wait first
 while (time-$time<3) {sleep 1}
 $time=time;
 
 my ($uid, $url, $qparam)=@_;
 $$qparam{'retmax'}=(scalar keys %$uid)+5;
 foreach my $ui (values (%$uid))
 {
  my $link=$url .  join "&", map {"$_=$$qparam{$_}"} keys %$qparam;
  $link .= "&id=$ui";
  $results.=LWP::Simple::get($link);
 }
 return $results;
}


sub tax_data {
 my $taxid=shift;
 my $link="http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=taxonomy&retmode=xml&id=$taxid";
 my $result=LWP::Simple::get($link);
 my @result=split /\n/, $result;
 my %result; # great, now we have a scalar, array and hash all called result!
 foreach (@result)
 {
  if (m#<Item Name="(.*?)" Type.*>(.*?)</Item>#)
  {
   $result{$1}=$2;
  }
 }
 # these are the things that we can get from the docsummary
 # Rank Division ScientificName CommonName TaxId NucNumber ProtNumber StructNumber GenNumber GeneNumber
 return ($result{'Rank'}, $result{'Division'}, $result{'ScientificName'}, $result{'CommonName'})
} 
  
