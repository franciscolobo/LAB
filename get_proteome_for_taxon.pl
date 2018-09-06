use strict;
use warnings;
use LWP::UserAgent;
use HTTP::Date;

#file containing genbank ids
my $id_file = $ARGV[0];

#chomp $top_node;

#store ids
my @ids;
my $i = 0;

open(IN, id_file);

while (my $line = <IN>) {
  chomp $line;
  $ids[$i] = $line;
  $i++;
}

close IN;

my $agent = LWP::UserAgent->new;

# Get a list of all reference proteomes of organisms below the given taxonomy node.
my $query_list = "http://www.uniprot.org/proteomes/?query=reference:no+taxonomy:$top_node&format=list";
my $response_list = $agent->get($query_list);
die 'Failed, got ' . $response_list->status_line .
  ' for ' . $response_list->request->uri . "\n" 
  unless $response_list->is_success;

# For each proteome, mirror its set of UniProt entries in compressed tabular format.
for my $proteome (split(/\n/, $response_list->content)) {
  my $file = $proteome . '.tab.gz';
  my $query_proteome = "http://www.uniprot.org/uniprot/?query=proteome:$proteome&columns=id,organism-id,proteome,go-id,database(KEGG)&format=tab&compress=yes";
  my $response_proteome = $agent->mirror($query_proteome, $file);

  if ($response_proteome->is_success) {
    my $results = $response_proteome->header('X-Total-Results');
    my $release = $response_proteome->header('X-UniProt-Release');
    my $date = sprintf("%4d-%02d-%02d", HTTP::Date::parse_date($response_proteome->header('Last-Modified')));
    print "File $file: downloaded $results entries of UniProt release $release ($date)\n";
  }
  elsif ($response_proteome->code == HTTP::Status::RC_NOT_MODIFIED) {
    print "File $file: up-to-date\n";
  }
  else {
    die 'Failed, got ' . $response_proteome->status_line .
      ' for ' . $response_proteome->request->uri . "\n";
  }
}
