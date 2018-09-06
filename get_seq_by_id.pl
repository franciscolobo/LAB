use strict;
use warnings;
use Bio::DB::GenBank;
use Bio::SeqIO;

my $usage = "perl get_seq ID out_dir\n\n";

if (!$ARGV[1]) {
  die ($usage);
}

my $id = $ARGV[0];
my $out_dir = $ARGV[1];
chomp $out_dir;

my $outfile = "$out_dir/$id.gb";

my $db_obj = Bio::DB::GenBank->new;

if (-e "$outfile") {
  print "$outfile already exists, skipping download\n";
} else {
  print ("Downloading $id\n");
  my $seq_obj = $db_obj->get_Seq_by_id($id);
  my $seqio_obj = Bio::SeqIO->new(-file => ">$outfile", -format => 'genbank' );
  $seqio_obj->write_seq($seq_obj);
}
