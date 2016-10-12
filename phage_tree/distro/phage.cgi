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

# phage script to get data from the databases.
# this will just allow people to get data.

# we will start with just getting all the genomes, and all the data from the genomes.

use strict;
use DBI;
use CGI::Carp qw(fatalsToBrowser);
use CGI qw/:standard/;
my $query = new CGI; # this is for the input data
my $output = new CGI; #this is for the output data
my $remaddr = $ENV{REMOTE_ADDR};
my $dbh = DBI->connect('DBI:mysql:phage') or croak "Can't connect to database\n";

# retrieve the parameters from the query. If there are none, then we will make the search table

my @param = $query->param;
unless (@param) {&blank_form()}    

if ($query->param('Submit') eq "BLAST") {&blast_form()}
if ($query->param('count') && $query->param('dbs')) {&get_by_count()}

# these are the current variables, note that I use a simple get_rows routine to get all
# the rows from the database, I can then change the rows as I like, and everything will be fine!
#PHAGE: organism, accession, locus, beginning, end, sequence
#PROTEIN: organism, start, stop, complement, gene, function, product, proteinid, dbxref, translation

unless (($query->param("genome") eq "ANY GENOME") || ($query->param("genome") eq "-----------")) {
	# we want some data from a genome
	# lets figure out what we want
	my $dbquery = 'SELECT count,organism,';
	my @columns = ('SEQ NUMBER', 'ORGANISM');
	if ($query->param('accession')) {$dbquery .= "accession,"; push (@columns, 'ACCESSION')}
	if ($query->param('locus')) {$dbquery .= "locus,"; push (@columns, "LOCUS\t")}
	if ($query->param('length')) {$dbquery .= 'beginning,end,'; push (@columns, ('BEGINNING', "END\t"))}
	if ($query->param('sequence')) {$dbquery .= 'sequence,'; push (@columns, 'SEQUENCE')}
	
	$dbquery =~ s/,$//;
	$dbquery .= " FROM phage WHERE ";
	
	if ($query->param("genome") eq "ALL GENOMES") {$dbquery .= "organism LIKE '%%'"} 
		else {$dbquery .= "organism LIKE '%".$query->param("genome")."%'"}
		
	unless ($query->param("keyword") eq "Keyword search") {
		my $rows = &get_rows('phage');
		$dbquery .= " AND (";
		my $keyword = " like '%".$query->param("keyword")."%' ||";
		$dbquery .= join ($keyword, @$rows).$keyword;
		$dbquery =~ s/ \|\|$/\)/;
		}
		
	my $outputtext ="<b>".$query->param('genome') . "</b><hr>";

	$outputtext .= &get_dna("$dbquery");
	&make_html("$outputtext");
	}
else {
	if ($query->param("keyword") eq "Keyword search") {&blank_form()}
	else {
		my $rows = &get_rows('phage');
		my $dbquery = "SELECT * FROM phage WHERE (";
		my $keyword = " like '%".$query->param('keyword')."%' || ";
		$dbquery .= join ($keyword, @$rows). $keyword;
		
		$dbquery =~ s/ \|\| $/\)/;
		$dbquery =~ s/ \|\| sequence like \'\%.*?\%\'//;
		
		my $outputtext = "<b>PHAGE GENOME MATCHES</b><p>";
		$outputtext .= &get_dna("$dbquery");
		
		
		# develop the query for the genome protein database
		$rows = &get_rows('protein');
		my $protquery = "SELECT * FROM protein WHERE (";
		$keyword = " like '%".$query->param('keyword')."%' || ";
		$protquery .= join ($keyword, @$rows). $keyword;
		$protquery =~ s/ \|\| $/\)/;
		$protquery =~ s/ \|\| translation like \'\%.*?\%\'//;
		
		#develop the query for the swiss protein database
		$rows = &get_rows('swiss');
		my $swissquery = "SELECT * FROM swiss WHERE (";
		$keyword = " like '%".$query->param('keyword')."%' || ";
		$swissquery .= join ($keyword, @$rows). $keyword;
		$swissquery =~ s/ \|\| $/\)/;
		$swissquery =~ s/ \|\| seq like \'\%.*?\%\'//;

		$outputtext .= &get_proteins("$protquery", "$swissquery", "2");

		
	&make_html("$outputtext");
	}
}



$dbh->disconnect;

sub get_all_orgs {
	my $exc = $dbh->prepare("SELECT organism FROM phage") or croak $dbh->errstr;
	$exc->execute or die $dbh->errstr;
	
	my @orgs = ('ANY GENOME', 'ALL GENOMES');
	while (my @retrieved = $exc->fetchrow_array) {push (@orgs, $retrieved[0])}
	return \@orgs;
}


sub get_rows {
	my $table=shift;
	my $exc = $dbh->prepare("show columns FROM $table") or croak $dbh->errstr;
	$exc->execute or die $dbh->errstr;
	
	my @rows;
	while (my @retrieved = $exc->fetchrow_array) {push (@rows, $retrieved[0])}
	return \@rows;
}

sub get_proteins {
	my $protsearch=shift;
	my $swisssearch=shift;
	my $tabs=shift;
	unless ($tabs) {$tabs="\t"} else {$tabs = "\t" x $tabs}

	
	my $temptext;
	if ($protsearch) {
	if ($query->param('protgenome')) {
		my $rows = &get_rows('protein');
		my $exc = $dbh->prepare("$protsearch") or croak $dbh->errstr;
		$exc->execute or die $dbh->errstr;
		
		$temptext .= "\n<b>PROTEINS FROM THE GENOME(S)</b>\n";
		while (my @retrieved = $exc->fetchrow_array) {
			foreach my $x (0 .. $#$rows) {
				if (uc($$rows[$x]) eq "TRANSLATION") {next unless ($query->param('protseq'))}
				next if (uc($$rows[$x]) eq "ORGANISM");
				next if (uc($$rows[$x]) eq "COUNT");
				next unless ($retrieved[$x]);
				if (uc($$rows[$x]) eq "PROTEINID") {
					$temptext .= $tabs.uc($$rows[$x])."\t<a href=\"http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Search&db=Protein&doptcmdl=GenPept&term=";
					$temptext .= $retrieved[$x]."\" target=\"_blank\">".$retrieved[$x]."</a>\n";
					next;
					}
				if ((uc($$rows[$x]) eq "DBXREF") && ($retrieved[$x] =~ /GI/)) {
					$temptext .= $tabs.uc($$rows[$x])."\t\t<a href=\"http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Search&db=Protein&doptcmdl=GenPept&term=";
					$retrieved[$x] =~ s/GI://; 
					$temptext .= $retrieved[$x]."\" target=\"_blank\">GI:".$retrieved[$x]."</a>\n";
					next;
					}
				if (uc($$rows[$x]) eq "TRANSLATION") {$retrieved[$x] =~ s/(.{40})/$1\n$tabs\t\t/g}
				if (uc($$rows[$x]) eq "COMPLEMENT") {chomp ($temptext); $temptext .= "\t(COMPLEMENT)\n"; next}
				if (length($$rows[$x]) < 8) {$temptext .= $tabs.uc($$rows[$x])."\t\t".$retrieved[$x]."\n"}
					else {$temptext .= $tabs.uc($$rows[$x])."\t".$retrieved[$x]."\n"}
			}
		$temptext .= "\n<hr>\n";
		}
	}
	}
	
	if ($swisssearch) {
	if ($query->param('protswiss')) {
		my $rows = &get_rows('swiss');
		my $exc = $dbh->prepare("$swisssearch") or croak $dbh->errstr;
		$exc->execute or die $dbh->errstr;
		
		$temptext .= "\n<b>SWISS PROT RESULTS</b>\n<hr>";
		while (my @retrieved = $exc->fetchrow_array) {
			foreach my $x (0 .. $#$rows) {
				next if (uc($$rows[$x]) eq "COUNT");
				# modify the output so that it fits nicely on the page
				if (uc($$rows[$x]) eq "SEQ") {$retrieved[$x] =~ s/(.{40})/$1\n$tabs\t\t/g}
					else {
					$retrieved[$x] =~ s/(.{50,60}.*?\s)/$1\n$tabs\t\t/g;
					if ($retrieved[$x] =~ /\-{10}/) {$retrieved[$x] =~ s/(\-{10,})/\n$tabs\t\t$1\n$tabs\t\t/g}
					$retrieved[$x] =~ s/\n\t+(.{0,20})\n/$1\n/g;
					$retrieved[$x] =~ s/^\s+//;
					}
				if (uc($$rows[$x]) eq "SEQ") {next unless ($query->param('protseq'))}
				next unless ($retrieved[$x]);
				if (uc($$rows[$x]) eq "PROTEINID") {
					$temptext .= $tabs.uc($$rows[$x])."\t<a href=\"http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Search&db=Protein&doptcmdl=GenPept&term=";
					$temptext .= $retrieved[$x]."\" target=\"_blank\">".$retrieved[$x]."</a>\n";
					next;
					}
				if ((uc($$rows[$x]) eq "DBXREF") && ($retrieved[$x] =~ /GI/)) {
					$temptext .= $tabs.uc($$rows[$x])."\t\t<a href=\"http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Search&db=Protein&doptcmdl=GenPept&term=";
					$retrieved[$x] =~ s/GI://; 
					$temptext .= $retrieved[$x]."\" target=\"_blank\">GI:".$retrieved[$x]."</a>\n";
					next;
					}
					
				if (uc($$rows[$x]) eq "COMPLEMENT") {chomp ($temptext); $temptext .= "\t(COMPLEMENT)\n"; next}
				if (length($$rows[$x]) < 8) {$temptext .= $tabs.uc($$rows[$x])."\t\t".$retrieved[$x]."\n"}
					else {$temptext .= $tabs.uc($$rows[$x])."\t".$retrieved[$x]."\n"}
			}
		$temptext .= "\n<hr>\n";
		}
	}
	}
	
	$temptext =~ s/<hr>\n$//;
	return $temptext;
}

sub get_dna {
	my $search=shift or die "can't get a search string\n";
	my $rows = &get_rows('phage');
	my $exc = $dbh->prepare("$search") or croak $dbh->errstr;
	$exc->execute or die $dbh->errstr;
	my $temptext;
	while (my @retrieved = $exc->fetchrow_array) {
		my $proteintext;
		foreach my $x (0 .. $#$rows) {
			if ($query->param('protein')) {
				my $search = "SELECT * from protein where organism = $retrieved[0]";
				$proteintext = &get_proteins("$search", "", "3")}
			if (uc($$rows[$x]) eq "SEQUENCE") {
				$retrieved[$x] =~ s/(.{60})/$1\n\t\t/g;
				$temptext .= $proteintext; undef $proteintext;
				}
			if (uc($$rows[$x]) =~ "ACCESSION") {
				$retrieved[$x] =~ s/GI\://;
				$temptext .= uc($$rows[$x])."\t<a href=\"http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Search&db=Nucleotide&doptcmdl=GenBank&term=";
				$temptext .= $retrieved[$x]."\" target=\"_blank\">GI:".$retrieved[$x]."</a>\n";
				next;
			}
			if (uc($$rows[$x]) =~ "LOCUS") {
				$temptext .= uc($$rows[$x])."\t\t<a href=\"http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Search&db=Nucleotide&doptcmdl=GenBank&term=";
				$temptext .= $retrieved[$x]."\" target=\"_blank\">".$retrieved[$x]."</a>\n";
				next;
			}
			$temptext .= uc($$rows[$x])."\t".$retrieved[$x]."\n";
			}
		if ($proteintext) {$temptext .= $proteintext; undef $proteintext}
		$temptext .= "<hr>"
		}
	return $temptext;
}
	
sub blank_form {
	my $orgs = &get_all_orgs();
	print $output -> header('text/html'),
                        start_html(     -title=>'Phage sequence databank',
                                                -BGCOLOR=>'#FFFFFF',
                                                -TEXT=>'#000000'),
        html((h1(big('Phage Sequence Databank'))),p,hr,p
	"Use this form to retrieve data from the phage sequence databank, or to perform a BLAST search (below) on the data",p,
	"Currently the phage sequence databank comes from two sources, the ",
	"<a href=\"http://www.ncbi.nlm.nih.gov:80/PMGifs/Genomes/phg.html\">Complete Phage Genomes</a>",
	" and the <a href=\"http://www.ebi.ac.uk/swissprot/\">SWISS PROT</a> protein database.",br,
	"Thanks to those sites for their data which I have extracted and cleaned up somewhat",p,
	"You can get all the data from a single phage, or search for keywords, etc.",p,hr,
	"To search the genome data:",p,
	startform(), popup_menu('genome', $orgs), p,hr,
	"Or you can do a keyword search. This will search all data except the actual protein or DNA sequences.",p,
	textfield('keyword','Keyword search',50),p,
	checkbox('protswiss', 'checked', 'ON', 'Search SWISS PROT phage'),
	checkbox('protgenome', 'checked', 'ON', 'Search genome data'),p,hr,
	"Retrieve this data: ",p,
	checkbox('accession','checked','ON','accesion number'), 
	checkbox('locus','checked','ON','locus'), 
	checkbox('length','checked','ON','sequence length'), 
	checkbox('sequence','','ON','DNA sequence'), 
	checkbox('protein','checked','ON','all proteins in the sequence'), 
	checkbox('protseq','','ON','sequences for the proteins'), p,
	hidden('count', ''), hidden ('dbs', ''),
	submit('Submit', 'SEARCH'), reset, p, hr, 
	h1("<b>BLAST SEARCH</b>"),p, "Note that the BLAST search is unaffected by the values above",p,
	radio_group('blastprogram', ['blastn', 'blastp', 'blastx', 'tblastn']),p,
	textarea('uploadsequence', '', 20, 80),p,
	submit('Submit', 'BLAST'), reset, p,
	p,end_form, small("Request came from $remaddr"));
	&exit_file('0');
	}


sub blast_form {

	my $blastprogram = $query->param('blastprogram');
	my @sequence = $query->param('uploadsequence');
	my $database;


#	if ($sequence =~ /\r/) {$sequence =~ s/\r/\n/g}
	unless ($sequence[0] =~ /^>/) {
		print $output -> header('text/html'),
			start_html(-title=>'ERROR', -BGCOLOR=>'#FFFFFF', -TEXT=>'#FF0000'),
		 html((h1(big('ERROR'))),p,hr,p
				 "The sequence not in FASTA format.",p,
				 "Please press the back button and try again by adding a line beginning > to the sequence",p
				 "If this is a problem, contact Rob Edwards (redwards\@utmem.edu)"),
		 end_html;
	 &exit_file('-1');
	 }
	
	if (($blastprogram eq "blastn") || ($blastprogram eq "tblastn") ||($blastprogram eq "tblastx")) {$database = "phage.nt.dbs"} else {$database = "phage.aa.dbs"}
	
	open (DNA, ">/tmp/temp$$") or die;

	#this obscure line sets the output for the text file to be immediate and returns all other output to as it was.
	select((select(DNA),$|=1)[0]);

	#The sequence needs to always be in fasta format
	print DNA "@sequence\n";
	close DNA;

	open(RESULTS, "nice -5 /usr/local/genome/bin/$blastprogram $database /tmp/temp$$  2>&1|"); 

	my $tempoutput;
	while (<RESULTS>) {$tempoutput .= $_}

	close (RESULTS);

	if ($tempoutput =~ /FATAL:  Could not read query sequence/) {
		 print $output -> header('text/html'),
				 start_html(	-title=>'ERROR',
							 -BGCOLOR=>'#FFFFFF',
							 -TEXT=>'#FF0000'),
		 html((h1(big('ERROR'))),p,hr,p
				 "There was an error with the sequence.",p,
				 "Please press the back button and try again",p,
				 "If this is a problem, contact Rob Edwards (redwards\@utmem.edu)"),
		 end_html;
	 &exit_file('-1');
	 }
	 
	 $tempoutput =~ /Sequences producing High-scoring Segment Pairs(.*?)>/s;
	 my $matches = $1;
	 my @hits;
	while ($matches =~ /\n\S+/) {
	 	$matches =~ s/\n(\S+)//;
		push (@hits, $1);
		}
	 
	print STDERR "MATCHES: ", join ("x", @hits), "\n\n";
	
	foreach my $hit (@hits) {
		my $resultstext; my $temp;
		if ($hit =~ /\d_\d/) {
			# it is from the protein database. Split it, and get the parts we need
			my ($count, $org) = split (/_/, $hit);
			my $exc = $dbh->prepare("SELECT gene,function,product,note from protein where count = $count") or croak $dbh->errstr;
			$exc->execute or die $dbh->errstr;
			my $done;
			while (my @retrieved = $exc->fetchrow_array) {
				foreach my $retrieved (@retrieved) {if ($retrieved) {$resultstext = $retrieved; $done=1; last}}
			}
			unless ($done) {$resultstext = "unknown"}
			my $earlyresults = $resultstext;
			
			$resultstext = "<a href=\"/cgi-bin/phage.cgi?&count=$count&dbs=protein&accession=1&locus=1&length=1&protein=1&protgenome=1&protseq=1\">". $resultstext."</a>";
			$exc = $dbh->prepare("SELECT organism from phage where count = $org") or croak $dbh->errstr;
			$exc->execute or die $dbh->errstr;
			while (my @retrieved = $exc->fetchrow_array) {
				$temp = $retrieved[0];
				$temp =~ s/\s+/\+/g;
				$resultstext .= " from <a href=\"/cgi-bin/phage.cgi?&count=$org&dbs=phage&accession=1&locus=1&length=1&protein=1&protgenome=1\">". $retrieved[0]."</a>";
				}
			$temp.=$earlyresults; # this is so it is the right length for the replacement later.
			}
		elsif ($hit =~ /__/) {
			# it is from swiss prot, so get the appropriate data from there
			my ($count, $ac) = split (/__/, $hit);
			my $exc = $dbh->prepare("SELECT gn,ac from swiss where count = $count") or croak $dbh->errstr;
			$exc->execute or die $dbh->errstr;
			my $done;
			while (my @retrieved = $exc->fetchrow_array) {
				foreach my $retrieved (@retrieved) {if ($retrieved) {$resultstext = $retrieved; $done=1; last}}
			}
			unless ($done) {$resultstext = "SWISS: unknown"}
			$resultstext =~ s/\.$//; 
			$resultstext =~ s/^\s+//;
			$resultstext = "SWISS: $resultstext";
			$temp = $resultstext;
			$resultstext = "<a href=\"/cgi-bin/phage.cgi?&count=$count&dbs=swiss&protswiss=1\">$resultstext</a>";
			}
			
		else {
			# it is from the phage (DNA) database
			
			my $exc = $dbh->prepare("SELECT organism from phage where count = $hit") or croak $dbh->errstr;
			$exc->execute or die $dbh->errstr;
			while (my @retrieved = $exc->fetchrow_array) {$resultstext = $retrieved[0]}
			$temp = $resultstext;
			$temp =~ s/\s+/\+/g;
			$resultstext = "<a href=\"/cgi-bin/phage.cgi?&count=$hit&dbs=phage&accession=1&locus=1&length=1&protein=1&protgenome=1\">". $resultstext."</a>";
			}
		$tempoutput =~ s/\n>$hit(\s+)/\n>$resultstext$1/g;
		my $diff =length($temp)-length($hit); # hack to remove the right number of spaces
		if ($diff <=0) {$diff = 1}
		$tempoutput =~ s/(Sequences producing High-scoring Segment Pairs.*?\n)$hit\s{0,$diff}/$1$resultstext/s;
		}

	 print $output -> header('text/html'),
				 start_html(	-title=>'Blast Results',
							 -BGCOLOR=>'#FFFFFF'),
				 html((h1("Blast Results")),p,hr,p
					 "<pre> $tempoutput </pre> \n",p),
				 end_html;


	#tidy up and exit
	unlink "/tmp/temp$$";
	&exit_file('0');
}


sub exit_file {
	my $exit= shift;
	$dbh->disconnect;
	exit($exit);
	}



sub make_html {
	my $outputtext=shift;
		print $output -> header('text/html'),
		start_html(-title=>'Phage sequence databank data', -BGCOLOR=>'#FFFFFF',  -TEXT=>'#000000'),
        	html((h1(big('Phage Sequence Databank DATA'))),p,hr,p
		"This is the data that you requested. All links will open in a new window. Thanks to the appropriate databases for the links.",p,hr,
		"<pre>$outputtext</pre><hr>",
		small("Request came from $remaddr"));
		&exit_file('0');
	}
	
sub get_by_count {
	my $count = $query->param('count');
	my $dbs = $query->param('dbs');
	my $outputtext;
	if ($dbs eq "phage") {$outputtext = &get_dna("SELECT * from phage where count = $count")}
	if ($dbs eq "swiss") {$outputtext = &get_proteins("", "SELECT * from swiss where count = $count", "1")}
	if ($dbs eq "protein") {$outputtext = &get_proteins("SELECT * from protein where count = $count", "2")}
print STDERR "For $count and $dbs got |$outputtext|\n";
	&make_html("$outputtext");
	}
		
