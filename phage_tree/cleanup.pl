$file = shift;
open (IN, $file) or die "Can't open $file\n";

while (<IN>) {
	next if (/\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*/);
	next if (/To Obtain Help Documentation:  send e-mail to /);
	next if (/with the word 'help' in the body of the mail message/);
	next if (/This is ENTREZ-Query server, that use daily updatable Entrez service/);
	next if (/Requested WWW query/);
	next if (/http:\/\/www.ncbi.nlm.nih.gov\/htbin-post\/Entrez\/query/);
	next if (/Content-type: text\/html/);
	next if (/Entrez Reports/);
	next if (/----------------/);
	print;
	}
