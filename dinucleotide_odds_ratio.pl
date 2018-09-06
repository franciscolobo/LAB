#!/usr/bin/perl

################################################################################
##                                                                            ##
## Copyright 2018 Universidade Federal de Minas Gerais                        ##
## Authors: Tarcisio Jose Domingos Coutinho/Francisco Pereira Lobo            ##
## this program is free software: you can redistribute it and/or modify       ##
## it under the terms of the GNU General Public License as published by the   ##
## Free Software Foundation, version 3 of the License.                        ##
##                                                                            ##
## extract_ORFs.pl is distributed in the hope that it will be useful,         ##
## but WITHOUT ANY WARRANTY; without even the implied warranty of             ##
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                       ##
## See the GNU General Public License for more details.                       ##
##                                                                            ##
## You should have received a copy of the GNU General Public License          ##
## along with rscu.pl (file: COPYING).                                        ##
##                                                                            ##
## If not, see <http://www.gnu.org/licenses/>.                                ##
##                                                                            ##
################################################################################

##takes as input a fasta file containing CDS data and prints
# genomic signature of individual sequences or for all sequences
# in this file

use strict;
use warnings;
use Bio::SeqIO;

#perl bin/dinucleotide_odds_ratio.pl $infile $out_dir_path $log_file_path $data{$key}{file_path} $strandedness 0 $data{$key}{group}

my $infile = $ARGV[0];

my $out_dir_path = $ARGV[1];

my $log_file_path = $ARGV[2];

my $flag = $ARGV[3];

my $mode = $ARGV[4];

my $label = $ARGV[5];

chomp $label;

if (($flag == 0) || ($flag == 1)) {

}
else {
  die ("Fourth argument ($flag) must be 0 (ss nucleic acid) or 1 (ds nucleic acid)\n");
}

if (($mode == 0) || ($mode == 1)) {
}
else {
  die ("Fifth argument ($mode) must be 0 (whole-genome) or 1 (individual sequences)\n");
}

my $file = "";
$file = get_fancy_name($infile);
my $outfile = "$out_dir_path$file\_dinucleotide_odds_ratio.txt";

if (-e $outfile) {
  exit(1);
} else {
  open (OUT, ">$outfile");
}

my %dinuc_count; #count data for each dinucleotide. Primary keys are dinucs (mode 0) or 
  #gene ids (mode 1)

my %nuc_count; #count data for each nucleotide. Primary keys as above;

my %stats; #final values

my $i = 0;

if ($mode == 0) {
  my ($tmp_hash_1, $tmp_hash_2) = create_dinuc_structured_data(0);
  %dinuc_count = %{$tmp_hash_1};
  %nuc_count = %{$tmp_hash_2};
} elsif ($mode == 1) {
  my ($tmp_hash_1, $tmp_hash_2) = create_dinuc_structured_data(1, "$infile");
  %dinuc_count = %{$tmp_hash_1};
  %nuc_count = %{$tmp_hash_2};
}

my $total_nts = 0; #total count of nts (mode 0)
my $total_dinucs = 0; #total count of dinucs (mode 0)

my $in = new Bio::SeqIO(-format => 'fasta',
                        -file => "$infile");
  
while (my $tmp_seq = $in->next_seq()) {
  my $sequence = $tmp_seq->seq();
  my ($tmp_1, $tmp_hash) = dinucleotide_freq($sequence);
  foreach my $key (keys %{$tmp_hash}) {
    $dinuc_count{$key} = $dinuc_count{$key} + $tmp_hash->{$key};
  }
  $total_dinucs = $total_dinucs + $tmp_1;
  $tmp_1 = 0;
  $tmp_hash = 0;
  ($tmp_1, $tmp_hash) = nuc_count($sequence);
  foreach my $key (keys %{$tmp_hash}) {
    $nuc_count{$key} = $nuc_count{$key} + $tmp_hash->{$key};
  }
  $total_nts = $total_nts + $tmp_1; 
}

if ($flag == 0) {
  print OUT ("GENOME\tLABEL\tDINUC\tCOUNT\tOBSERVED\tEXPECTED\tODDS\n");
  foreach my $key (keys %dinuc_count) {
    my $expected = expected_frequency(0, 0, $key, \%nuc_count, $total_nts); #returns expected frequency for dinucletide $key
    my $observed = $dinuc_count{$key} / $total_dinucs;
    my $odds = $observed / $expected;
    print OUT ("$file\t$label\t$key\t$dinuc_count{$key}\t$observed\t$expected\t$odds\n");
  #  print OUTFIN ("$genome\t$key\t$dinuc_count{$key}\t$observed\t$expected\t$odds\n");
  }
} 
elsif ($flag ==1) {
  print OUT ("GENOME\tLABEL\tDINUC\tCOUNT_DINUC\tREVERSE_DINUC\tCOUNT_REVERSE_DINUC\tOBSERVED\tEXPECTED\tODDS\n");
  foreach my $key (keys %dinuc_count) {
    my $expected = expected_frequency(1, 0, $key, \%nuc_count, $total_nts); #returns expected frequency for dinucletide $key
    my $reverse_key = revcom($key); #complement reverse
    my $observed = 2 * (($dinuc_count{$key} / $total_dinucs) + ($dinuc_count{$reverse_key}/$total_dinucs));
    my $odds = $observed / $expected;
    print OUT ("$file\t$label\t$key\t$dinuc_count{$key}\t$reverse_key\t$dinuc_count{$reverse_key}\t$observed\t$expected\t$odds\n");
  }
}

close OUT;

sub revcom {
  my $tmp_seq = $_[0];
  my $revcomp = reverse($tmp_seq);
  $revcomp =~ tr/ACGTacgt/TGCAtgca/;
  return $revcomp;
}

sub expected_frequency {
  my $tmp_flag = shift @_;
  my $tmp_mode = shift @_;
  if ($tmp_mode == 0) {
    my $dinuc = shift @_;
    my $tmp_hash = shift @_;
    my $tmp_total = shift @_;
    if ($tmp_flag == 0) { #ss dna
      my @aux = split (//, $dinuc);
      my $first_nt = shift @aux;
      my $second_nt = shift @aux;
      my %tmp = %{$tmp_hash};
      my $first_nt_count = $tmp{$first_nt};
      my $second_nt_count = $tmp{$second_nt};
#      print "$first_nt\t$first_nt_count\t$second_nt\t$second_nt_count\t$tmp_total\n\n";
#      my $a = <STDIN>;
      my $first_nt_freq = $first_nt_count / $tmp_total;
      my $second_nt_freq = $second_nt_count / $tmp_total;
      my $expected_freq = $first_nt_freq * $second_nt_freq;
      return ($expected_freq);
    }
    if ($tmp_flag == 1) {#ds dna
      my @aux = split (//, $dinuc);
      my $first_nt = shift @aux;
      my $second_nt = shift @aux;
      my %tmp = %{$tmp_hash};
      my $first_nt_count = $tmp{$first_nt};
      my $second_nt_count = $tmp{$second_nt};
      my $first_nt_count_rev = $tmp{revcom($first_nt)};
      my $second_nt_count_rev = $tmp{revcom($second_nt)};
      my $first_nt_freq = $first_nt_count / $tmp_total;
      my $second_nt_freq = $second_nt_count / $tmp_total;
      my $first_nt_rev_freq = $first_nt_count_rev / $tmp_total;
      my $second_nt_rev_freq = $second_nt_count_rev / $tmp_total;
      my $expected_freq = (($first_nt_freq + $second_nt_freq) * ($first_nt_rev_freq + $second_nt_rev_freq));
      return ($expected_freq);
    }
  }
  if ($tmp_mode == 1) {
    my $dinuc = shift @_;
    my $seq_id = shift @_;
    my $tmp_hash = shift @_;
    my $tmp_total = shift @_;
    if ($tmp_flag == 0) { #ss dna
      my @aux = split (//, $dinuc);
      my $first_nt = shift @aux;
      my $second_nt = shift @aux;
      my %tmp = %{$tmp_hash};
      my %tmp_total = %{$tmp_total};
      my $first_nt_count = $tmp{$seq_id}{$first_nt};
      my $second_nt_count = $tmp{$seq_id}{$second_nt};
      my $first_nt_freq = $first_nt_count / $tmp_total{$seq_id};
      my $second_nt_freq = $second_nt_count / $tmp_total{$seq_id};
      my $expected_freq = $first_nt_freq * $second_nt_freq;
      return ($expected_freq);
    }
    if ($tmp_flag == 1) {#ds dna
      my @aux = split (//, $dinuc);
      my $first_nt = shift @aux;
      my $second_nt = shift @aux;
      my %tmp = %{$tmp_hash};
      my %tmp_total = %{$tmp_total};
      my $first_nt_count = $tmp{$seq_id}{$first_nt};
      my $second_nt_count = $tmp{$seq_id}{$second_nt};
      my $first_nt_count_rev = $tmp{$seq_id}{revcom($first_nt)};
      my $second_nt_count_rev = $tmp{$seq_id}{revcom($second_nt)};
      my $first_nt_freq = $first_nt_count / $tmp_total{$seq_id};
      my $second_nt_freq = $second_nt_count / $tmp_total{$seq_id};
      my $first_nt_rev_freq = $first_nt_count_rev / $tmp_total{$seq_id};
      my $second_nt_rev_freq = $second_nt_count_rev / $tmp_total{$seq_id};
      my $expected_freq = (($first_nt_freq + $second_nt_freq) * ($first_nt_rev_freq + $second_nt_rev_freq));
      return ($expected_freq);
    }
  }
}

sub nuc_count {
  my $sequence = $_[0];
  my %tmp_hash;
  my @aux;
  @aux = split (//, $sequence);
  my $total = 0;
  my ($a, $c, $t, $g) = 0;
  foreach my $element (@aux) {
    if ($element eq "A") {
      $a++;
      $total++;
    }
    elsif ($element eq "C") {
      $c++;
      $total++;
    }
    elsif ($element eq "T") {
      $t++;
      $total++;
    }
    elsif ($element eq "G") {
      $g++;
      $total++;
    }
    else {
      die ($!);
      $total++;
    }
  }
  $tmp_hash{"A"} = $a;
  $tmp_hash{"C"} = $c;
  $tmp_hash{"T"} = $t;
  $tmp_hash{"G"} = $g;
  return ($total, \%tmp_hash);
}

sub dinucleotide_freq {
  my @DNA = ("A","C","T","G");
  my %dinuc;
  my $total;
  foreach my $letter (@DNA) {
    foreach my $letter_2 (@DNA) {
      my $actual = join ("", $letter, $letter_2);
      $dinuc{$actual} = 0;
    }
  }
  my $sequence = $_[0];
  my @aux;
  @aux = split (//, $sequence);
  for (my $i = 0; $i < $#aux; $i++) {
    my $actual = join ("", $aux[$i],$aux[$i+1]);
    $dinuc{$actual}++;
    $total++;
  }
  return ($total, \%dinuc);
}


sub create_dinuc_structured_data {
  my $tmp_mode = $_[0];
  my @letters = ("A", "C", "T", "G");
  my (%tmp_dinuc_count, %tmp_nuc_count);
  if ($tmp_mode == 0) {
    for (my $i = 0; $i <= $#letters; $i++) {
      $tmp_nuc_count{$letters[$i]} = 0;
      for (my $i2 = 0; $i2 <= $#letters; $i2++) {
        my $dinuc = join ("", $letters[$i], $letters[$i2]);
        $tmp_dinuc_count{$dinuc} = 0;
      }
    }
    return (\%tmp_dinuc_count, \%tmp_nuc_count);
  } elsif ($tmp_mode ==1) {
    my $tmp_infile = $_[1] ;
    my $in = new Bio::SeqIO(-format => 'fasta',
                            -file => $tmp_infile);
    my @all_ids;
    my $i = 0;
    while (my $tmp_seq = $in->next_seq()) {
      $all_ids[$i] = $tmp_seq->display_id();
      $i++;
    }
    foreach my $id (@all_ids) {
      for (my $i = 0; $i <= $#letters; $i++) {

        $tmp_nuc_count{$id}{$letters[$i]} = 0;
        for (my $i2 = 0; $i2 <= $#letters; $i2++) {
          my $dinuc = join ("", $letters[$i], $letters[$i2]);
          $tmp_dinuc_count{$id}{$dinuc} = 0;
        }
      }
    }
    return (\%tmp_dinuc_count, \%tmp_nuc_count);
  }
}

sub get_fancy_name {
  my $tmp = $_[0];
  my @aux = split(/\//, $tmp);
  $tmp = pop @aux;
  @aux = "";
  @aux = split(/\./, $tmp);
  $tmp = shift @aux;
  $tmp =~ s/_nt$//g;
  return $tmp;
}
