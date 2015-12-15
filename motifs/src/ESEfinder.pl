#!/bin/perl

use strict;
use Bio::Tools::Analysis::DNA::ESEfinder;
use Bio::SeqIO;

my $in = shift @ARGV;
my $out = shift @ARGV;
open(my $o, '>', $out);
print $o "motif\tsr_protein\tscore\n";

my $seqio = Bio::SeqIO->new(-file => "$in", '-format' => 'Fasta', '-alphabet' => 'dna');
while(my $seq = $seqio->next_seq) {
  my $string = $seq->seq;
  print "$string\n";
  my $pseq = Bio::PrimarySeq->new(
       -id => 'tmp',
       -seq=>"$string",
       -alphabet=>'dna');
  my $ese_finder = Bio::Tools::Analysis::DNA::ESEfinder->new(-seq => $pseq);
  $ese_finder->run();
  die "Could not get a result"
    unless $ese_finder->status =~ /^COMPLETED/;
 
  foreach my $feat ( $ese_finder->result('Bio::SeqFeatureI') ) {
    my $motif = ($feat->get_tag_values('motif'))[0];
    my $protein = ($feat->get_tag_values('SR_protein'))[0];
    my $score = ($feat->get_tag_values('score'))[0];
    print $o "$motif\t$protein\t$score\n";
  }
}
close $o;

