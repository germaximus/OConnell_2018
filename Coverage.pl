#!/usr/bin/perl
# calculates coverage for every gene in the bowtie output.
# Also requires mRNA_100uniq.fasta containing list of genes (bowtie index)
# Usage: Coverage.pl ARGV[0]
# where ARGV[0] - is the uniq.bwt output of a bowtie-1

my %genes;

#initialize index hash. Due to memory requirements I have to run it on a server
open (INDEX, "<mRNA_100uniq.fasta") || die "no index file found";
while ($line = <INDEX>)
  {  chomp($line);
     $name = substr($line,1);
     $line = <INDEX>; chomp($line);
     $length = length($line);
     for(my $i=0; $i < $length; $i++)    {   $genes{$name}[$i] = 0;    }
  }
close (INDEX);


open (INPUT, "<$ARGV[0]"); open (OUT, ">coverage"); 
while ($line = <INPUT>)
  {
     @line = split /\t/,$line;
     
     for(my $i = $line[3]; $i < ($line[3] + length($line[4])); $i++)
       {
         $genes{$line[2]}[$i]++;
       }
  }
foreach $name (sort keys %genes)
   {
      print OUT "$name\t";
      print OUT join("\t",@{$genes{$name}})."\n";
   }
close(INPUT); close(OUT);
