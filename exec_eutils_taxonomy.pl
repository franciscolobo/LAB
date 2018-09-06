use strict;
use warnings;

my $infile = $ARGV[0];

chomp $infile;

open(IN,"<$infile");

my @queries;

my $i = 0;

while (my $line = <IN>) {
  chomp $line;
  $queries[$i] = $line;
  $i++;
}

close IN;

foreach my $query (@queries) {
  my $tmp = $query;
  $query =~ s/\s+/ /g;
  $tmp =~ s/\s+/_/g;
  my $species_filename = "$tmp.xml";
  if (-e $species_filename) {
    print "$species_filename already downloaded\n";
    next;
  } else {
    print "./esearch -db taxonomy -query \"$query\[orgn\]\" | ./efetch -format txt > $species_filename\n";
#  system "./esearch -db taxonomy -query \"$query\[taxid\]\" | ./efetch -format txt > $tmp.txt\n";
    system "./esearch -db taxonomy -query \"$query\[orgn\]\" | ./efetch -format xml > $species_filename\n";
  }
#  system "./esearch -db nucleotide -query \"$query\[orgn\] AND biomol_mrna\[PROP\] AND complete CDS\[title\] NOT mitochondrial\[title]\" | ./efetch -format gb > $tmp.gb\n";
}
