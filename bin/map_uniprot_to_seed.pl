#__perl

use strict;
use SAPserver;
use Data::Dumper;
my $sap=new SAPserver;
use Getopt::Std;
my %opts;
getopts('f:n:', \%opts);

die <<EOF unless ($opts{f});
$0
-f filename to get
-n number of proteins to parse at once. The default is 1,000 proteins


This is a program to retrieve essential seed data (role, subsystems, and subsystem classifications) about a set of proteins. To map into the seed we map based on the MD5 sum of the uppercase protein sequence (just the protein sequence alone, and without any identifiers or stop symbol). We can then move from those ids to subsystems and roles.

This code uses the seed servers, and you can find more information about them at http://servers.theseed.org/ You will ned to install the seed servers before starting, but that should install all the commands you need.

Basically, we make remote calls that retrieve data, and the calls are designed so that we can retrieve a bolus of data and compute on it. In this case, you can set how many pieces of data to send/receive with the -n option. I suggest a relatively large number, like 1000, since the delay is typically in retrieving the data and not in sending the data over the internet.

This will print the results out, and so you may need to redirect the output.

EOF


unless (defined $opts{n}) {$opts{n} = 1000}

# get the subsytem hierarchy. We only need this once, so getting it here will let us use it later
my $ssHash = $sap->all_subsystems(-usable => 1);

# counters and raw data
my $ids; my $count=0;
# open the file and iterate through every line
# this allows us to open gzipped or zipped files
if ($opts{f} =~ /gz$/) {
	open(IN, "gunzip -c $opts{f}|") || die "Can't open a pipe to gunzip $opts{f}";
} elsif ($opts{f} =~ /zip$/) {
	open(IN, "unzip -p $opts{f}|") || die "can't open a pipe to unzip $opts{f}";
}       else {
	open(IN, $opts{f}) || die "can't open $opts{f}";
}

while (<IN>) {
	chomp;
	# this is the input data
	my ($id, $ac, $go, $md5, $trans)= split /\t/;
	# store it and make a call with lots of data
	$ids->{$md5}=[$id, $ac, $go, $md5, $trans];
	$count++;
	if ($count == $opts{n}) {
		# we need to make the call and print the information out
		# don't forget to reset our counters.
		&get_data($ids);
		$count=0;
		$ids={};
	}
}
close IN;


sub get_data {
	my $ids = shift;

	# convert the md5 sums to fig ids
	# all of the subsequent calls require fig ids
	my $idHash = $sap->proteins_to_fids({-prots=>[keys %$ids]});
	my $figids;
	foreach my $id (keys %$idHash) {
		map {$figids->{$_}=1} @{$idHash->{$id}};
	}

	# get the subsystems for each proteins
	my $subsysHash = $sap->ids_to_subsystems({-ids=>[keys %$figids], -usable=>1});
	
	# now print out everything
	# we need to iterate through the md5s and the fid ids to merge everything together
	foreach my $md5 (keys %$ids) {
		foreach my $fid (@{$idHash->{$md5}}) {
			# get the functional role and subsystem for this protein
			# recall that this is a one:many relationship
			foreach my $tple (@{$subsysHash->{$fid}}) {
				my ($role ,$ss)=@$tple;
				# in case the classification isn't defined, we'll make it null
				if (!defined $ssHash->{$ss}) {$ssHash->{$ss}=['', ['', '']]}
				# print everything out
				print join("\t", @{$ids->{$md5}}, $fid, $role, $ss, @{$ssHash->{$ss}->[1]}), "\n";
			}
		}
	}
}



