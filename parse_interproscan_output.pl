use strict;
use warnings;

if (!$ARGV[2]) {
  print_help();
} 

my $infile = $ARGV[0];

my $outfile = $infile."_parsed.txt";

my $column_number = $ARGV[1];

my $feature_name = $ARGV[2];

chomp $feature_name;

open(IN, "<$infile") || die($!);

print OUT ("Feature\t$feature_name\n");

my $i = 0;

while (my $line = <IN>) {
  chomp $line;   
  my @aux = split(/\t/, $line);
  my $actual_feature = $aux[$column_number-1];
  if (defined $actual_feature) {
    print OUT ("$i\t$actual_feature\n");
    $i++;
  }
}

close IN;

sub print_help {
  die("Use this program like: perl parse_interproscan_tabular.pl <path to interproscan output file> <column to summarize> <Feature name to be used as header>\n");
}

