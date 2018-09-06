use strict;
use warnings;
use Bio::SeqIO;

my $flag = 0; #to compute for codons, change to 1

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

my $tmp;

$nrow = $#seqs;

if ($flag == 1) {
  for (my $i2 = 0; $i2 <= $#seqs; $i2++) {
    my $seq = $seqs[$i2];
    my @aux = split (//, $seq);
    $ncol = (($#aux+1)/3)-1;
    for (my $i = 0; $i <= $ncol; $i++) {
      my $codon = "";
      my $first = shift @aux;
      my $second = shift @aux;
      my $third = shift @aux;
      $codon = $first.$second.$third;
#      print "$codon\n";
#      my $a = <STDIN>;
      $values[$i2][$i] = $codon; #[seq_index][char_index]
    }
    $tmp = $i;
  }
} else {
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
}

#print "$ncol\t$nrow\t$tmp\n\n";
#my $a = <STDIN>;

my @results;

my $total = 0;

for (my $i = 0; $i <= $ncol; $i++) {
  my %count;
  for (my $i2 = 0; $i2 <= $nrow; $i2++) {
    if (defined $count{$values[$i2][$i]}) {
      $count{$values[$i2][$i]}++;    
    }
    else {
      $count{$values[$i2][$i]} = 1;
    }
  }
  my $sum = 0;
  my $total = $nrow+1;
  foreach my $key (keys %count) {
    my $frequency = $count{$key}/($total);
    my $log_result = log2(1/$frequency);
    $sum = $sum + ($frequency*$log_result);
  }
  print "$sum\n";
}


sub log2 {
  my $n = shift;
  return log($n)/log(2);
}
