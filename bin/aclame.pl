#!/usr/bin/perl -w

# Web whack aclame

use strict;
use LWP::UserAgent;
use HTTP::Request::Common qw(POST GET);
use LWP::Simple;

# http://aclame.ulb.ac.be/perl/Aclame/Tools/export_clusters.cgi
# 


my %dx=(
    'type' => "cluster",
    'fmt' => "fasta",
    'vers' => 0.2,
    'cat' => "all",
    'class' => "family",
    'submit' => "FASTA"
);

#    'id' => "cluster:all:0",

for (my $i=0; $i<=427; $i++)
{
    my $url='http://aclame.ulb.ac.be/perl/Aclame/Tools/export_clusters.cgi?';
    map {$url .= "$_=$dx{$_}&"} keys %dx;
    $url .= "id=cluster:all:$i";
    print STDERR "$url\n";
    print get($url), "\n";
}
exit(0);


my $ua=LWP::UserAgent->new;
$ua->agent("Rob_Edwards/0.1 ");

my $req = GET 'http://aclame.ulb.ac.be/perl/Aclame/Tools/export_clusters.cgi',
[
    'type' => "cluster",
    'fmt' => "fasta",
    'vers' => 0.2,
    'cat' => "all",
    'id' => "cluster:all:0",
    'class' => "family",
    'submit' => "FASTA"
];


print $ua->request($req)->as_string;

