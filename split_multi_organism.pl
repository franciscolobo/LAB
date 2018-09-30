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


use strict;
use warnings;
use Bio::SeqIO;

my $usage = "\nperl split_multi_organism_gb.pl <genbank file containing more than one species>\n\n";

if (!$ARGV[0]) {
  die ($usage);
}

my $infile = $ARGV[0];
my $in = new Bio::SeqIO(-format => 'genbank',
                        -file => $infile);


my @tmp = split (/\//, $infile);

my $trash = pop @tmp;

my $outdir = join ("/", @tmp);


my $seq;
my %organisms;

while($seq = $in->next_seq()){
  for my $feat_object ($seq->get_SeqFeatures) {
#    print "primary tag: ", $feat_object->primary_tag, "\n";
    for my $tag ($feat_object->get_all_tags) {
      if ($tag eq "organism") {
        for my $specie ($feat_object->get_tag_values($tag)) {
#          print ("$specie\n");
          if (defined $organisms{$specie}) {
#            print ("Opa!\n");
  	  }
	  else {
	    $organisms{$specie} = 1;
	  }
        }
      }
      for my $value ($feat_object->get_tag_values($tag)) {
#         print "    value: ", $value, "\n";
      }
    }
  }
}

foreach my $key (keys %organisms) {
  my $outfile = $key;
  $outfile =~ s/\s+/_/g;
  $outfile =~ s/-/_/g;
  $outfile =~ s/\\/_/g;
  $outfile =~ s/\//_/g;
  $outfile =~ s/\(/_/g;
  $outfile =~ s/\)/_/g;
  print ("$outfile\n");
  my $in = new Bio::SeqIO(-format => 'genbank',
                          -file => $infile);
  my $seqout = Bio::SeqIO->new(-file   => ">$outdir/$outfile.gb",
                               -format => 'genbank');

  my $seq;
  while ($seq = $in->next_seq()){
    for my $feat_object ($seq->get_SeqFeatures) {
#    print "primary tag: ", $feat_object->primary_tag, "\n";
      for my $tag ($feat_object->get_all_tags) {
        if ($tag eq "organism") {
          for my $specie ($feat_object->get_tag_values($tag)) {
            print ("\t->$specie\n");
            if ($key eq $specie) {
              print ("\t\t!!!$key\n");
              $seqout->write_seq($seq);
            }
          }
        }
      }
    }
  }
} 
