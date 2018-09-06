use strict;
use warnings;
use Bio::SeqIO;
use Parallel::ForkManager;


my $pm = new Parallel::ForkManager(4);

my $indir = shift @ARGV;

my @files;

opendir(DIR, "$indir");

@files = readdir(DIR);

closedir(DIR);

foreach my $file (@files) {
#  my $pid = $pm->start and next; 
  my $file_path = $indir.$file;
  next if ($file !~ /.gb_nt.fasta$/);
#  $pm->finish if ($file !~ /.gb_nt.fasta$/);
  my $tmpname = $indir.$file.".single_gene.fasta";
  print "$tmpname\n";
#  my $a = <STDIN>;
  if (-e $tmpname) {
    print "Here!\n";
    next;
  }
  my $seqio_object = Bio::SeqIO->new(-file => $file_path); 
  my $superseq; #
  while (my $seq_object = $seqio_object->next_seq) {
    my $seq = $seq_object->seq();
    $superseq = $superseq.$seq;
  }
  my @aux = split(/\./, $file);
  my $name = shift @aux;
  my $outfile = $indir.$file.".single_gene.fasta";
  open (OUT, ">$outfile");
  print OUT ">$name\n$superseq\n";
  close OUT;
  $pm->finish; # do the exit in the child process
}
$pm->wait_all_children;