#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;

chomp $ARGV[0] if defined;

chomp $ARGV[1] if defined;

if ((!defined $ARGV[1])||(!-d $ARGV[0])||(!-d $ARGV[1])||($ARGV[0] eq "-h")||(($ARGV[0] eq "--help"))) { #checking if users are asking for help
  help();
  exit;
}

my $dir_path = shift;

my $output_dir = shift;

opendir(DIR,$dir_path);

my @files = readdir (DIR);

foreach my $file (@files) {
  if (($file eq ".") || ($file eq "..")) {
    next;
  }
  if ($file !~ /.gb$/) {
    next;
  }
  print ("$file\n");
#  my $a = <STDIN>;
  my $path_2_file = "$dir_path$file";
  my $in = new Bio::SeqIO(-format => 'genbank',
                           -file => $path_2_file);
  my $seq;

  my $name = $file;
  my $i = 1;

  open (OUT, ">", "$output_dir/$file\_proteins.fasta");
  open (OUT_2, ">", "$output_dir/$file\_nt.fasta");
  open (LOG, ">", "$output_dir/$file\_IDs.log");

  while($seq = $in->next_seq()){
    my $sequence = $seq->seq();
    my $start;
    my $end;
    my $organism;
    my $source;
    my $cds;
    my $flag = 0;
    my $desc = $seq->desc();
    my $strain = "";
#    print "$desc\n";
    for my $feat_object ($seq->get_all_SeqFeatures) {
      print "primary tag: ", $feat_object->primary_tag, "\n";
      for my $tag ($feat_object->get_all_tags) {
        for my $value ($feat_object->get_tag_values($tag)) {
          if ($tag eq "pseudo") { #remove pseudogenes
#            print ("$tag\t$value\n");
            $flag = 1;
#            next;
          }
          if ($tag eq "strain") {
            print "$tag = $value\n";
            $strain = $value;
#            print "Opa!\n";
#            my $a = <STDIN>;
          }
        }
      }
     if ( $feat_object->primary_tag eq 'CDS' ) {
       my $cds_object = $feat_object->spliced_seq;
       $cds = $cds_object->seq;
       if ($cds !~ /^ATG/) {
#       print ("$name\t$organism\_$i\t$cds\n");
        }
#       print "CDS sequence is ",$cds_obj->seq,"\n";
     }
       if ($feat_object->primary_tag eq "source") {
         $source = $file;
         $source =~s/\.gb.*$//g;
         $source =~ s/\s+/_/g;
#        print "primary tag: ", $feat_object->primary_tag, "\n";
         for my $tag ($feat_object->get_all_tags) {
#          print ("$tag\n");
          if ($tag eq "organism") {
            for my $value ($feat_object->get_tag_values($tag)) {
              $organism = $value;
            }
          }
        }
      $organism =~ s/\s+/_/g;
      $organism =~ s/\-/_/g;
      }
      if ($feat_object->primary_tag eq "CDS") {
        if ($flag == 1) {
          $flag = 0;
          next;
        }
#       print "primary tag: ", $feat_object->primary_tag, "\n";
        print OUT (">$source||$desc||$strain||$i\n");
        print OUT_2 (">$source||$desc||$strain||$i\n$cds\n");
        for my $tag ($feat_object->get_all_tags) {
          if ($tag eq "protein_id") {
            for my $value ($feat_object->get_tag_values($tag)) {
              print LOG ("$source||$desc||$strain||$i\t$value\n");
              $i++;
            }
          }
         if ($tag eq "translation") {
            for my $value ($feat_object->get_tag_values($tag)) {
              print OUT ("$value\n");
            }
          }
        }
#      $organism =~ s/\s+/_/g;
      }
    }
  }
}

sub help {
  print "\nperl extract_aa_nt_from_gb.pl <directory with genbank files> <directory to write output>\n\n";
  print "The script will read gb files, extract CDS features and translations and write CDS and protein data in two fasta files\n\n";
  print "A third file, linking the common temporary ID used in both files to the old protein/CDS ids is also created\n\n";
  exit();
}
