use strict;
use warnings;
use Bio::SeqIO;

my $infile = $ARGV[0];

my $outdir = $ARGV[1];

my $label = $ARGV[2];

chomp $label;

my $name = fancy_name($infile);

my $outfile_path = $outdir.$name."_nuc_per_codon_pos.txt";

compute_data_nucpos($infile, $outfile_path);

sub compute_data_nucpos {
  my $infile = shift @_;
  my $outfile = shift @_;
  open(OUT,">$outfile");
  print OUT ("GENOME\tLABEL\tNUCLEOTIDE\tPOSITION\tVALUE\n");
  my $in = new Bio::SeqIO(-format => 'fasta',
                          -file => "$infile");
  
  my $hash = "";
  my %tmp_data;
  $hash = create_nuc_structure(\%tmp_data);
  %tmp_data = %{$hash};
  while (my $seq = $in->next_seq()) {
    my $tmp_seq = $seq->seq();
    my @aux = '';
    @aux = split (//, $tmp_seq);
    my $i = $#aux;
    while ($i >= 0) {
      my $actual_codon = join ('', $aux[0], $aux[1], $aux[2]);
      $tmp_data{$aux[0]}{1}++;
      $tmp_data{total}{1}++;
      $tmp_data{$aux[1]}{2}++;
      $tmp_data{$aux[2]}{3}++;
      shift @aux;
      shift @aux;
      shift @aux;
      $i = $i - 3;
    }
  }
  foreach my $key (keys %tmp_data) {
    my @sec_keys = ("1", "2", "3");
    for (my $i = 0; $i <= $#sec_keys; $i++) {
      my $sec_key = $sec_keys[$i];
      my $freq = $tmp_data{$key}{$sec_key}/$tmp_data{total}{1};
      next if ($key eq "total");
      print OUT "$name\t$label\t$key\t$sec_key\t$freq\n";
    }
  }
  close OUT;
}

sub create_nuc_structure {
  my $tmp_data = shift @_;
  my %tmp_hash = %{$tmp_data};
  my @nucs = ("A", "C", "T", "G");
  my @poss = ("1", "2", "3");
  foreach my $nuc (@nucs) {
    foreach my $pos (@poss) {
      $tmp_hash{$nuc}{$pos} = 0;
    }
  }
  foreach my $pos (@poss) {
    $tmp_hash{total}{$pos} = 0;
  }
  return \%tmp_hash;
}

sub fancy_name {
  my $tmp = $_[0];
  my @aux = split(/\t/, $tmp);
  my @aux2 = split(/\//, $aux[0]);
  my $file_name = pop @aux2;
  $file_name =~ s/.fasta$//;
  $file_name =~ s/_nt$//;
  return $file_name;
}

