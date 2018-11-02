#!/usr/bin/perl
# modifies *.gff3 annotation file by dropping off "Curated Genomic" pseudogenes and Gnomon records. gff file has to have its header.

open (INPUT, "<GRCh38.p12.custom.gff"); 
for($i=0; $i < 8; $i++) {$line = <INPUT>; print $line;}
while ($line = <INPUT>)
  {
     @fields = split /\t/, $line;
     if($fields[1] ne 'Gnomon' && $fields[1] ne 'Curated Genomic')  {  print $line;  }
  }
close(INPUT); 