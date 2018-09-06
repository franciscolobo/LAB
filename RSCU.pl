#!/usr/bin/perl

################################################################################
##                                                                            ##
## Copyright 2018 Universidade Federal de Minas Gerais                        ##
## Author: Tarcisio Jose Domingos Coutinho/Francisco Pereira Lobo             ##
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

# Takes as input a genbank file and calculates RSCU for all genes
# after filtering putative error sources.

use strict;
use warnings;
use Bio::SeqIO;
use Bio::Tools::CodonTable;

#if (!$ARGV[1]) {
#  print ("use this program like: perl <genbank file> <0 or 1 to indicate whole-genome or individual gene output format>\n");
#  die;
#}

my $infile = $ARGV[0];

my $outfile_path = $ARGV[1];

my $logfile_path = $ARGV[2];

my $codon_table = $ARGV[3];

my $label = $ARGV[4];

my $mode = $ARGV[5]; #0 for whole-genome, 1 for individual values for each gene

chomp $mode;

my @dna = ("A", "C", "T", "G");
my %number_of_codons2amino; #degeneracy number for each aminoacid
my %codon2amino; #translation
my %all_codons; #all codons
my %codon2values; #count data for each analysis
my $e = 0; #index
my $total_codons = 0; #totao count of codons
if ($mode == 0) {
  my  ($tmp_hash_1, $tmp_hash_2, $tmp_hash_3, $tmp_hash_4) = generate_codon_structured_data(0, $codon_table); #will populate data structures for mode 0
  %all_codons = %{$tmp_hash_1};
  %codon2amino = %{$tmp_hash_2};
  %codon2values = %{$tmp_hash_3};
  %number_of_codons2amino = %{$tmp_hash_4};
} elsif ($mode == 1) {
  my ($tmp_hash_1, $tmp_hash_2, $tmp_hash_3, $tmp_hash_4) = generate_codon_structured_data(1, $codon_table, $infile); #will populate data structures for mode 1
  %all_codons = %{$tmp_hash_1};
  %codon2amino = %{$tmp_hash_2};
  %codon2values = %{$tmp_hash_3};
  %number_of_codons2amino = %{$tmp_hash_4};
} else {
  die(); #mode must be equal to 0 or 1
}
#now I generated three hash strucutures and one array structure:
# %number_of_codons2amino - has as key a one-letter aminoacid code and returns the number of codons for it
# %codon2amino - has as key a given codon and returns the one letter code for it
# %codon2values - this hash was initialized here and will contain as key one codon and will return the number of times this codon was observed in a given sequence
# @all_codons and array wich contains all codons;
my $in = new Bio::SeqIO(-format => 'fasta',
                        -file => "$infile");
my %total_values_of_amino; #number of times an amino acid is coded.
my @all_codons = keys %all_codons;

if ($mode == 0) {
  for (my $i = 0; $i <= $#all_codons; $i++) { #initializing
      $total_values_of_amino{$codon2amino{$all_codons[$i]}} = 0;
  }
} elsif ($mode == 1) {
  while (my $seq = $in->next_seq()){
    for (my $i = 0; $i <= $#all_codons; $i++) { #initializing
      my $acc = $seq->display_id();
      $total_values_of_amino{$acc}{$codon2amino{$all_codons[$i]}} = 0;
    }
  }
}
my @all_ids; #used to store all genes ids if analysis mode eq 1

my $i2 = 0; #used as index for @all_ids if analysis mode eq 1

my %total_codons_per_gene; #used in analysis mode eq 1

$in = new Bio::SeqIO(-format => 'fasta',
                     -file => "$infile");
if ($mode == 0) {
  my $genome = $infile;
  $genome =~ s/\.gb_nt\.fasta$//g;
  $genome =~ s/^tmp\///g;
  while (my $seq = $in->next_seq()){
    my $sequence = $seq->seq();
    my $acc = $seq->display_id();
#      print "$acc\n";
    my @aux = '';
    my $total_codons = 0;
    @aux = split (//, $sequence);
#      print "$sequence\n";
    my $i = $#aux;
    while ($i >= 0) {
      $total_codons = $total_codons + 1;
      my $actual_codon = join ('', $aux[0], $aux[1], $aux[2]);
#        print $actual_codon."\n";
#        my $a = <STDIN>;
      shift @aux;
      shift @aux;
      shift @aux;
      $codon2values{$actual_codon} = $codon2values{$actual_codon} + 1;
      $total_values_of_amino{$codon2amino{$actual_codon}} = $total_values_of_amino{$codon2amino{$actual_codon}} + 1;
      $i = $i - 3;
    }
  }
} elsif ($mode == 1) {
  while (my $seq = $in->next_seq()){
    my $sequence = $seq->seq();
    my $acc = $seq->display_id();
    my @aux = '';
    my $total_codons = 0;
    @aux = split (//, $sequence);
    my $i = $#aux;
    while ($i >= 0) {
      $total_codons = $total_codons + 1;
      my $actual_codon = join ('', $aux[0], $aux[1], $aux[2]);
      shift @aux;
      shift @aux;
      shift @aux;
#     print ("*$acc\t*$actual_codon*\t*$codon2values{$acc}{$actual_codon}*\t*$total_values_of_amino{$acc}{$codon2amino{$actual_codon}}*\t*$number_of_codons2amino{$codon2amino{$actual_codon}}*\n");
#        my $a = <STDIN>;
      $codon2values{$acc}{$actual_codon} = $codon2values{$acc}{$actual_codon} + 1;
      $total_values_of_amino{$acc}{$codon2amino{$actual_codon}} = $total_values_of_amino{$acc}{$codon2amino{$actual_codon}} + 1;
      $all_ids[$i2] = $acc;
      $i = $i - 3;
      $i2++;
    }
    $total_codons_per_gene{$acc} = $total_codons;
    $total_codons = 0;
  }
} else {
    die ("Analysis mode $ARGV[0] must be equal 0 or 1\n");
}

#getting a fancy name
my $genome = get_fancy_name($infile);
my $outfile = "$outfile_path$genome\_rel_syn_cod_usa.txt";
open(OUT, ">$outfile");

if ($mode == 0) {
  print OUT ("GENOME\tLABEL\tCODON\tAA\tCODON_COUNT\tAA_COUNT\tAA_DEGENERACY\tEXPECTED\tRSCU\n");
  foreach my $key (%codon2values) {
    if (exists $codon2values{$key}) {
      my $expected;
      if (exists $number_of_codons2amino{$codon2amino{$key}}) {
        if ($number_of_codons2amino{$codon2amino{$key}} == 0) {
          $expected = 0;
        }
        else {
          $expected = ($total_values_of_amino{$codon2amino{$key}}/$number_of_codons2amino{$codon2amino{$key}});
        }
      }
      else {
        $expected = 0;
      }
      my $codon_usage;
      if ($expected == 0) {
        $codon_usage = 0;
      }
      else {
        $codon_usage = ($codon2values{$key}/$expected);
      }
      print OUT ("$genome\t$label\t$key\t$codon2amino{$key}\t$codon2values{$key}\t$total_values_of_amino{$codon2amino{$key}}\t$number_of_codons2amino{$codon2amino{$key}}\t$expected\t$codon_usage\n");
    }
  }
} elsif ($mode == 1) {
  print OUT ("GENOME\tLABEL\tGENE_ID\tCODON\tAA\tCODON_COUNT\tAA_COUNT\tAA_DEGENERACY\tEXPECTED\tRSCU\n");
  foreach my $key (@all_ids) {
   foreach my $codon (@all_codons) {
      if (exists $codon2values{$key}{$codon}) {
        my $expected;
        if (exists $number_of_codons2amino{$codon2amino{$codon}}) {
          if ($number_of_codons2amino{$codon2amino{$codon}} == 0) {
            $expected = 0;
          }
          else {
            $expected = ($total_values_of_amino{$key}{$codon2amino{$codon}}/$number_of_codons2amino{$codon2amino{$codon}});
          }
        }
        else {
          $expected = 0;
        }
        my $codon_usage;
        if ($expected == 0) {
          $codon_usage = 0;
        }
        else {
          $codon_usage = ($codon2values{$key}{$codon}/$expected);
        }
      print OUT ("$genome\t$label\t$key\t$codon\t$codon2amino{$codon}\t$codon2values{$key}{$codon}\t$total_values_of_amino{$key}{$codon2amino{$codon}}\t$number_of_codons2amino{$codon2amino{$codon}}\t$expected\t$codon_usage\n");
      }
    }
  }
}

sub generate_codon_structured_data {
  my (%tmp_all_codons, %tmp_codon2amino, %tmp_codon2values, %tmp_number_of_codons2amino);
  my $tmp_mode = $_[0]; #analysis mode
  my $codon_table = $_[1];
  my @tmp_ids; #will store ids if mode eq 1;
  my $myCodonTable = Bio::Tools::CodonTable->new(-id => $codon_table);
  my $i = 0;
  if ($tmp_mode == 1) {
    my $tmp_infile = $_[2];
    my $in = new Bio::SeqIO(-format => 'fasta',
                            -file => $tmp_infile);

    while (my $seq = $in->next_seq()){
      $tmp_ids[$i] = $seq->display_id();
      $i++;
    }
  }
  my @dna = ("A", "C", "T", "G");
    for (my $i = 0; $i <= $#dna; $i++) {
    my $first_letter = $dna[$i];
    for (my $i_2 = 0; $i_2 <= $#dna; $i_2++) {
      my $second_letter = $dna[$i_2];
      for (my $i_3 = 0; $i_3 <= $#dna; $i_3++) {
        my $third_letter = $dna[$i_3];
        my $codon = join ('', $first_letter, $second_letter , $third_letter);
        $tmp_all_codons{$codon} = 1;
        $tmp_codon2amino{$codon} = $myCodonTable->translate($codon);
        if ($tmp_mode == 0) {
          $tmp_codon2values{$codon} = 0;
        }
        if ($tmp_mode == 1) {
          foreach my $tmp_id (@tmp_ids) {
#            print "$tmp_id\t$codon\n";
#            my $a = <STDIN>;
            $tmp_codon2values{$tmp_id}{$codon} = 0
          }
        }
        if (defined $tmp_number_of_codons2amino{$tmp_codon2amino{$codon}}) {
          $tmp_number_of_codons2amino{$tmp_codon2amino{$codon}} = $tmp_number_of_codons2amino{$tmp_codon2amino{$codon}} + 1;
        }
        else {
          $tmp_number_of_codons2amino{$tmp_codon2amino{$codon}} = 1;
        }
      }
    }
  }
  return (\%tmp_all_codons, \%tmp_codon2amino, \%tmp_codon2values, \%tmp_number_of_codons2amino);
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

