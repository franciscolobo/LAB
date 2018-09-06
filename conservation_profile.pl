use strict;
use warnings;
use Bio::SeqIO;

my $io = Bio::SeqIO->new(-file   => "$ARGV[0]",
                               -format => "fasta" );


my @seqs;
my @ids;
my $i = 0;
while (my $seqobj = $io->next_seq) { 
  my $id = $seqobj->display_id();
  my $seq = $seqobj->seq;
  $seqs[$i] = $seq;
  $ids[$i] = $id;
  $i++;
#  print ("$id\t$seq\n");
}

my @values;

$i = 0;

my $nrow = 0;
my $ncol = 0;

$nrow = $#seqs;

for (my $i2 = 0; $i2 <= $#seqs; $i2++) {
  my $seq = $seqs[$i2];
  my @aux = split (//, $seq);
  $ncol = $#aux;
  foreach my $element(@aux) {
    $values[$i2][$i] = $element;
    $i++;
  }
  $i = 0;
}


my @results;

for (my $i = 0; $i <= $ncol; $i++) {
  my %count;
  for (my $i2 = 0; $i2 <= $nrow; $i2++) {
#    print ("$i\t$i2\t$values[$i2][$i]\n");
    if (defined $count{$values[$i2][$i]}) {
      $count{$values[$i2][$i]}++;    
    }
    else {
      $count{$values[$i2][$i]} = 1;
    }
  }
  my $highest = 0;
  my $letter = "";
  foreach my $key (keys %count) {
    if ($count{$key} > $highest) {
      $highest = $count{$key};
      $letter = $key;
    }
  }
#  if ($letter eq "-") {
#    print ("$i\t$letter\t0\n");
#    next;
#  }
  print ("$i\t$letter\t$highest\t");
  my @aux;
  my $i = 0;
  foreach my $key (sort {$count{$a} cmp $count{$b}} (keys %count)) {
    $aux[$i] = "$key:$count{$key}";
    $i++;
  }
  my $tmp = join(";", @aux);
  print "$tmp\n";
#   my $a = <STDIN>;
}
