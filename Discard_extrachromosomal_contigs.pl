#!/usr/bin/perl

my %chromosomes =('NC_000001.11',0,
                  'NC_000002.12',0,
                  'NC_000003.12',0,
                  'NC_000004.12',0,
                  'NC_000005.10',0,
                  'NC_000006.12',0,
                  'NC_000007.14',0,
                  'NC_000008.11',0,
                  'NC_000009.12',0,
                  'NC_000010.11',0,
                  'NC_000011.10',0,
                  'NC_000012.12',0,
                  'NC_000013.11',0,
                  'NC_000014.9',0,
                  'NC_000015.10',0,
                  'NC_000016.10',0,
                  'NC_000017.11',0,
                  'NC_000018.10',0,
                  'NC_000019.10',0,
                  'NC_000020.11',0,
                  'NC_000021.9',0,
                  'NC_000022.11',0,
                  'NC_000023.11',0,
                  'NC_000024.10',0 );

open (INPUT, '<', '../NCBI_GRCh38.p12_originals/GRCh38.p12.fna') or die "Can't open file";


while ($line = <INPUT>)  {
     @line = split('\s+', $line);
     if(substr($line[0],0,1) eq '>' && exists($chromosomes{substr($line[0],1)})) {
           print $line;
           while ($line = <INPUT>) {
                  if (substr($line,0,1) ne '>') { print $line;   }
                  else {last;}                
           }
           redo; 
     }
    
}
close(INPUT);






