#use strict;
#use warnings;
#use warnings;
##Program that realizes a global alignment of two given sequences and prints as result a table, containing in each cell a score value and a pointer indicating the correct way of alignment. Based on The Needleman-Wunsch algorithm. Rum this program like: Use this program like: perl alinhamento_global_blosum.pl <first_sequence> <second_sequence> <flag indicating the usage of BLOSUM45 (1), BLOSUM50 (2) or BLOSUM62 (3)> <value for gap>

##Defining some The BLOSUM Tables and other sutff that will be essential. Those matrices were obtained from the web.

# The BLOSUM45 amino acid substitution matrix as a hash array.

@blosum45 =
    #  A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
  ( [  5,-2,-1,-2,-1,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-2,-2, 0],  # A
    [ -2, 7, 0,-1,-3, 1, 0,-2, 0,-3,-2, 3,-1,-2,-2,-1,-1,-2,-1,-2],  # R
    [ -1, 0, 6, 2,-2, 0, 0, 0, 1,-2,-3, 0,-2,-2,-2, 1, 0,-4,-2,-3],  # N
    [ -2,-1, 2, 7,-3, 0, 2,-1, 0,-4,-3, 0,-3,-4,-1, 0,-1,-4,-2,-3],  # D
    [ -1,-3,-2,-3,12,-3,-3,-3,-3,-3,-2,-3,-2,-2,-4,-1,-1,-5,-3,-1],  # C
    [ -1, 1, 0, 0,-3, 6, 2,-2, 1,-2,-2, 1, 0,-4,-1, 0,-1,-2,-1,-3],  # Q
    [ -1, 0, 0, 2,-3, 2, 6,-2, 0,-3,-2, 1,-2,-3, 0, 0,-1,-3,-2,-3],  # E
    [  0,-2, 0,-1,-3,-2,-2, 7,-2,-4,-3,-2,-2,-3,-2, 0,-2,-2,-3,-3],  # G
    [ -2, 0, 1, 0,-3, 1, 0,-2,10,-3,-2,-1, 0,-2,-2,-1,-2,-3, 2,-3],  # H
    [ -1,-3,-2,-4,-3,-2,-3,-4,-3, 5, 2,-3, 2, 0,-2,-2,-1,-2, 0, 3],  # I
    [ -1,-2,-3,-3,-2,-2,-2,-3,-2, 2, 5,-3, 2, 1,-3,-3,-1,-2, 0, 1],  # L
    [ -1, 3, 0, 0,-3, 1, 1,-2,-1,-3,-3, 5,-1,-3,-1,-1,-1,-2,-1,-2],  # K
    [ -1,-1,-2,-3,-2, 0,-2,-2, 0, 2, 2,-1, 6, 0,-2,-2,-1,-2, 0, 1],  # M
    [ -2,-2,-2,-4,-2,-4,-3,-3,-2, 0, 1,-3, 0, 8,-3,-2,-1, 1, 3, 0],  # F
    [ -1,-2,-2,-1,-4,-1, 0,-2,-2,-2,-3,-1,-2,-3, 9,-1,-1,-3,-3,-3],  # P
    [  1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-3,-1,-2,-2,-1, 4, 2,-4,-2,-1],  # S
    [  0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-1,-1, 2, 5,-3,-1, 0],  # T
    [ -2,-2,-4,-4,-5,-2,-3,-2,-3,-2,-2,-2,-2, 1,-3,-4,-3,15, 3,-3],  # W
    [ -2,-1,-2,-2,-3,-1,-2,-3, 2, 0, 0,-1, 0, 3,-3,-2,-1, 3, 8,-1],  # Y
    [  0,-2,-3,-3,-1,-3,-3,-3,-3, 3, 1,-2, 1, 0,-3,-1, 0,-3,-1, 5]   # V
    #  A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
    );

# The BLOSUM50 amino acid substitution matrix as a hash array.

@blosum50 =
    #  A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
  ( [  5,-2,-1,-2,-1,-1,-1, 0,-2,-1,-2,-1,-1,-3,-1, 1, 0,-3,-2, 0],  # A
    [ -2, 7,-1,-2,-4, 1, 0,-3, 0,-4,-3, 3,-2,-3,-3,-1,-1,-3,-1,-3],  # R
    [ -1,-1, 7, 2,-2, 0, 0, 0, 1,-3,-4, 0,-2,-4,-2, 1, 0,-4,-2,-3],  # N
    [ -2,-2, 2, 8,-4, 0, 2,-1,-1,-4,-4,-1,-4,-5,-1, 0,-1,-5,-3,-4],  # D
    [ -1,-4,-2,-4,13,-3,-3,-3,-3,-2,-2,-3,-2,-2,-4,-1,-1,-5,-3,-1],  # C
    [ -1, 1, 0, 0,-3, 7, 2,-2, 1,-3,-2, 2, 0,-4,-1, 0,-1,-1,-1,-3],  # Q
    [ -1, 0, 0, 2,-3, 2, 6,-3, 0,-4,-3, 1,-2,-3,-1,-1,-1,-3,-2,-3],  # E
    [  0,-3, 0,-1,-3,-2,-3, 8,-2,-4,-4,-2,-3,-4,-2, 0,-2,-3,-3,-4],  # G
    [ -2, 0, 1,-1,-3, 1, 0,-2,10,-4,-3, 0,-1,-1,-2,-1,-2,-3, 2,-4],  # H
    [ -1,-4,-3,-4,-2,-3,-4,-4,-4, 5, 2,-3, 2, 0,-3,-3,-1,-3,-1, 4],  # I
    [ -2,-3,-4,-4,-2,-2,-3,-4,-3, 2, 5,-3, 3, 1,-4,-3,-1,-2,-1, 1],  # L
    [ -1, 3, 0,-1,-3, 2, 1,-2, 0,-3,-3, 6,-2,-4,-1, 0,-1,-3,-2,-3],  # K
    [ -1,-2,-2,-4,-2, 0,-2,-3,-1, 2, 3,-2, 7, 0,-3,-2,-1,-1, 0, 1],  # M
    [ -3,-3,-4,-5,-2,-4,-3,-4,-1, 0, 1,-4, 0, 8,-4,-3,-2, 1, 4,-1],  # F
    [ -1,-3,-2,-1,-4,-1,-1,-2,-2,-3,-4,-1,-3,-4,10,-1,-1,-4,-3,-3],  # P
    [  1,-1, 1, 0,-1, 0,-1, 0,-1,-3,-3, 0,-2,-3,-1, 5, 2,-4,-2,-2],  # S
    [  0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 2, 5,-3,-2, 0],  # T
    [ -3,-3,-4,-5,-5,-1,-3,-3,-3,-3,-2,-3,-1, 1,-4,-4,-3,15, 2,-3],  # W
    [ -2,-1,-2,-3,-3,-1,-2,-3, 2,-1,-1,-2, 0, 4,-3,-2,-2, 2, 8,-1],  # Y
    [  0,-3,-3,-4,-1,-3,-3,-4,-4, 4, 1,-3, 1,-1,-3,-2, 0,-3,-1, 5]   # V
    #  A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
  );

# The BLOSUM62 amino acid substitution matrix as a hash array.

@blosum62 =
    #  A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
  ( [  4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0],  # A
    [ -1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3],  # R
    [ -2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3],  # N
    [ -2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3],  # D
    [  0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1],  # C
    [ -1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2],  # Q
    [ -1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2],  # E
    [  0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3],  # G
    [ -2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3],  # H
    [ -1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3],  # I
    [ -1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1],  # L
    [ -1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2],  # K
    [ -1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1],  # M
    [ -2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1],  # F
    [ -1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2],  # P
    [  1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2],  # S
    [  0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0],  # T
    [ -3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3],  # W
    [ -2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1],  # Y
    [  0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4]   # V
    #  A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
    );

##Starting the program and checking if everythingh is OK

# The two first arguments of the program are the two sequences to be aligned. The third indicates which BLOSUM matrice will be used and the fourth  is the gap value, respectively.

if (!$ARGV[3]) {
  die ("Use this program like: $0 <first_sequence> <second_sequence> <flag indicating the usage of BLOSUM45 (1), BLOSUM50 (2) or BLOSUM62 (3)> <value for gap>\n");
}

$gap = $ARGV[3]; ## Value for GAP;

##Defining which matrix will be used.

if (($ARGV[2]) == 1) {
  @matrix = @blosum45;
}
elsif (($ARGV[2]) == 2) {
  @matrix = @blosum50;
}
elsif (($ARGV[2]) == 3) {
  @matrix = @blosum62;
}
else {
  print ("\$ARGV[2] must have a value of 0, 1 or 2 to indicate the usage of BLOSUM 45, 50 or 62 respectivelly\n");
  die;
}

# Creating a array whith the amino acids.

@aas  = split(//, "ARNDCQEGHILKMFPSTWYV");

##Creating a hashmap %score_of_alignment that links each pair of aminoacids to it's alignment score. Like: A|R (this is the key) has -2 (this is the value) in the hashmap

for ($i=0; $i<20; $i++) {
  for ($j=0; $j<20; $j++) {
    $score_of_alignment{"$aas[$i]|$aas[$j]"} = $matrix[$i][$j]; #This is iscanning the matrix choosed by $ARGV[2] and passing it to th hashmap.
  }
}

$seq_1 = $ARGV[0]; #first sequence
$seq_2 = $ARGV[1]; #second sequence

##Translating everithing to an unique format of letters, removing whitespaces and checking if every letter is a valid aminoacid code;

$seq_1 =~ tr/[a-z]/[A-Z]/;
$seq_2 =~ tr/[a-z]/[A-Z]/;

$seq_1 =~ s/\s//g;
$seq_2 =~ s/\s//g;

if (($seq_1 =~ /\d/) || ($seq_1 =~ /[BJOUXZ]/)) {
  print ("The first sequence has non-aminoacid characters\n");
  die;
}

if (($seq_2 =~ /\d/) || ($seq_2 =~ /[BJOUXZ]/)) {
  print ("The second sequence has non-aminoacid characters\n");
  die;
}

##Splitting each sequence in it's letters.

@letter_seq_1 = split (//, $seq_1);
@letter_seq_2 = split (//, $seq_2);

$i = 0; #pointer to horizontal position;
$j = 0; #pointer to vertical position;

##staring to fill the first row and collum of matrix @score with the initial values of gaps;

$length_seq_1 = length($seq_1);
$length_seq_2 = length($seq_2);

$score[$i][$j] = 0;
$pointer[$i][$j] = "null";

## filling the first line of the atual vector

for ($i = 1; $i <= $length_seq_1; $i++) {
  $pointer[$i][0] = "-";
}

for ($i = 1; $i <= $length_seq_2; $i++) {
  $atual[0] = ($gap * $i);
  $pointer[0][$i] = "|";
}

#for ($a = 0; $a <= $length_seq_2; $a++) {
#  for ($b = 0; $b <= $length_seq_1; $b++) {
#    print ("$pointer[$b][$a]\t");
#  }
#  print ("\n");
#}

##filling the rest of the matrix; The two "for" will fill all the matrix with the correct values.

for ( $j = 1; $j <= $length_seq_2; $j++) {
  $letter_2 = $letter_seq_2[$j-1];
    @previous = @atual; ##the @atual will be cleaned and their values will be passed to @previous
    @atual = '';
    $atual[0] = $gap*$j; #filling the first cell with the gap value;
    for ( $i = 1; $i <= $length_seq_1; $i++) {
      $letter_1 = $letter_seq_1[$i-1];
      $diagonal_score = '';
      $left_score = '';
      $up_score = '';
      $diagonal_score = $previous[$i-1] + $score_of_alignment{"$letter_1|$letter_2"}; ##This will fill the diagonal score with the alignment value.
      $left_gap_score = $atual[$i-1] + $gap; #Same as digonal score
      $up_gap_score = $previous[$i] + $gap;  #Same as diagonal score
      if ($diagonal_score >= $up_gap_score) {
        if ($diagonal_score >= $left_gap_score) {
          $atual[$i] = $diagonal_score;
	  $pointer[$i][$j] = "\\";
#	  print ("$pointer[$i][$j]||$letter_1||$letter_2\t");
        }
        else {
          $atual[$i] = $left_gap_score;
	  $pointer[$i][$j] = "-";
#	  print ("$pointer[$i][$j]||$letter_1||$letter_2\t");
        }
     }
     elsif ($up_gap_score > $left_gap_score) {
       $atual[$i] = $up_gap_score;
       $pointer[$i][$j] = "|";
#       print ("$pointer[$i][$j]||$letter_1||$letter_2\t");
     }
     else {
       $atual[$i] = $left_gap_score;
       $pointer[$i][$j] = "-";
#       print ("$pointer[$i][$j]||$letter_1||$letter_2\t");

     }
     # The "if" and "else" structure above is checking the greatest value of a given comparisson of two letters, i.e., it's checking whta score (and it's pointer) is the gratest one and picking it as the chosen. In a case where the same value is found, the diagonal one is the chosen.
  }
#  print ("\n");
}

##in this point 

for ($a = 0; $a <= $length_seq_2; $a++) {
  for ($b = 0; $b <= $length_seq_1; $b++) {
    print ("\t$pointer[$b][$a]");
  }
  print ("\n");
}

$b = $length_seq_1;
$a = $length_seq_2;

if ($a <= $b) {
  $max_size = $b-1;
}
else {
  $max_size = $a-1;
}

until (($a == 0) && ($b ==0)) {
    $atual_pointer = $pointer[$b][$a];
    if ($atual_pointer eq "\\") {
      $seq_1_final[$max_size] = $letter_seq_1[$b-1];
      $seq_2_final[$max_size] = $letter_seq_2[$a-1];
#      print ("$a\t$b\t$atual_pointer\n");
#      $temp = join ("", @seq_1_final);
#      print ("1\t$temp\n");
#      $temp = join ("", @seq_2_final);
#      print ("2\t$temp\n");
      $a--;
      $b--;
      $max_size--;
    }
    elsif ($atual_pointer eq "-"){
      $seq_1_final[$max_size] = $letter_seq_1[$b-1];
      $seq_2_final[$max_size] = "-";
#      $temp = join ("", @seq_1_final);
#      print ("1\t$temp\n");
#      $temp = join ("", @seq_2_final);
#      print ("2\t$temp\n");
#      print ("$a\t$b\t$atual_pointer\n");
      $b--;
      $max_size --;
    }
    elsif ($atual_pointer eq "|") {
      $seq_1_final[$max_size] = "-";
      $seq_2_final[$max_size] = $letter_seq_2[$a-1];
#      $temp = join ("", @seq_1_final);
#      print ("1\t$temp\n");
#      $temp = join ("", @seq_2_final);
#      print ("2\t$temp\n");
#      print ("$a\t$b\t$atual_pointer\n");
      $a--;
      $max_size --;
    }
    else {
    }
#  $au = <STDIN>;
#  print ("\n\n");
}

$seq_1_string_final = join ("", @seq_1_final);
$seq_2_string_final = join ("", @seq_2_final);

print ("$seq_1_string_final\n");
print ("$seq_2_string_final\n");
##now we have the @atual vetor full and it's last value is th score of alignment.

#variable $#atual points to the end of vector. Build in function of perl.
#if you wish to print the sequences, please remove the "#" character of the next two lines;

$score = $atual[$#atual];
print ("The score value of alignment is $score\n");