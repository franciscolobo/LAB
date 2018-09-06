#!/usr/bin/perl

################################################################################
##                                                                            ##
## Copyright 2015 Embrapa Informatica Agropecuaria                            ##
## Author: Francisco Pereira Lobo                                             ##
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

my $dir_path = shift;

my $output_dir = shift;

my $len_cutoff = shift;

my $flags = shift;

if (! defined $len_cutoff) {
  $len_cutoff = 100;
}

if (! defined $flags) {
  $flags = "all";
}

opendir(DIR,$dir_path);

my @files = readdir(DIR);

foreach my $file (@files) {
  if (($file eq ".") || ($file eq "..")) {
    next;
  }
  if ($file !~ /.gb$/) {
    next;
  }
  print ("$file\n");
  my $path_2_file = "$dir_path$file";
  my $in = new Bio::SeqIO(-format => 'genbank',
                           -file => $path_2_file);
  my $seq;

  my $name = $file;
  my $i = 1;

  open (OUT, ">", "$output_dir/$file\_nt.fasta");
#  open (LOG, ">", "$output_dir/$file\_proteins.log");

  while($seq = $in->next_seq()){
    my $sequence = $seq->seq();
    my $start;
    my $end;
    my $organism;
    my $source;
    my $cds;
    my $flag = 0;
    for my $feat_object ($seq->get_all_SeqFeatures) {
#      print "primary tag: ", $feat_object->primary_tag, "\n";
      for my $tag ($feat_object->get_all_tags) {
        for my $value ($feat_object->get_tag_values($tag)) {
          if ($tag eq "pseudo") { #remove pseudogenes
            $flag = 1;
            next;
          }
        }
      if ( $feat_object->primary_tag eq 'CDS' ) {
        my $cds_object = $feat_object->spliced_seq;
        $cds = $cds_object->seq;
      }
      if ($feat_object->primary_tag eq "source") {
        $source = $file;
	$source =~s/\.gb.*$//g;
        $source =~ s/\s+/_/g;
        for my $tag ($feat_object->get_all_tags) {
          if ($tag eq "organism") {
            for my $value ($feat_object->get_tag_values($tag)) {
              $organism = $value;
            }
          }
        }
      $organism =~ s/\s+/_/g;
      $organism =~ s/\-/_/g;
      }
    }
      if ($feat_object->primary_tag eq "CDS") {
        next if check_cds($cds); #next if cds contains any problem
        print OUT (">$source\_$i");
        for my $tag ($feat_object->get_all_tags) {
          if ($tag eq "protein_id") {
            for my $value ($feat_object->get_tag_values($tag)) {
              print OUT ("||protein_id:$value");
              $i++;
            }
          }
          if ($tag eq "organism") {
            for my $value ($feat_object->get_tag_values($tag)) {
              print OUT ("||organism:$value");
            }
          }
          if ($tag eq "locus_tag") {
            for my $value ($feat_object->get_tag_values($tag)) {
              print OUT ("||locus_tag:$value");
            }
          }
          if ($tag eq "gene") {
            for my $value ($feat_object->get_tag_values($tag)) {
              print OUT ("||gene:$value");
            }
          }
        }
        print OUT ("\n$cds\n");
      }
    }
    $i++;
  }
}

sub check_cds {
  my $tmp_seq = $_[0];
  my $start_codon = substr($tmp_seq, 0, 3);
  my $stop_codon = substr($tmp_seq, -3);
  if ($start_codon !~ /ATG|GTG/i) {
    print "Not valid start codon\t$tmp_seq\n";
    return (1);
  }
  if ($stop_codon !~ /TAA|TAG|TGA/i) {
    print "Not valid stop codon\t$tmp_seq\n";
    return (1);
  }
  if (length($tmp_seq) % 3 == 0) {
  } else {
    print "length not multiple of three\t$tmp_seq\n";
    return (1);
  }
  if ($tmp_seq =~ /[^ACGT]/i) {
    print "non-standard nucleotides\t$tmp_seq\n";
    return (1);
  }
  if (length($tmp_seq) < $len_cutoff) {
    print "length smaller than cutoff\n";
    return (1);
  }
  return (0);
}

