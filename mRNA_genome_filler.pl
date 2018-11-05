#!/usr/bin/perl
##fills 5'UTRs and 3'UTRs shorter than 100 with nucleotides inferred from .gff annotation and chromosome sequences.
# temp3 file is taken as a result of mRNA_extractor.pl
open (INPUT, "<temp3"); open (OUT, ">mRNA_100"); open (NO_MATCH, ">no_Match");

# Read mRNA genomic coordinates from the GFF annotation
my %annotation;
open (ANNOTATION, "<GRCh38.p12.Refseq.coding.gff");
while (my $line = <ANNOTATION>) {
   my @line = split /\t/, $line;
   if ( $line[2] eq 'mRNA'  && $line[8] =~ m/Genbank:([^;,\.]+)/) {
       $gene_name = $1;
       $annotation{$gene_name} = [ $line[0], $line[3], $line[4], $line[6] ];   }
} close(ANNOTATION);
print "extracting mRNA coordinates from GFF annotation completed\n";
#---------------------------------------------------------------

# Build a hash with chromosome coordinates
my %genome;
open (GENOME, "<GRCh38.p12.custom.fna");
my $seq = "";
my $chrom_name = "";
while (my $line = <GENOME>) {
    chomp($line);
    if ($line =~ /^>/) {
        if (length($seq) > 0) { $genome{$chrom_name} = $seq; $seq = ""; }
        @line = split /\s/, $line;
        $chrom_name = substr($line[0], 1); }   
    else { $seq = $seq.$line}
}
$genome{$chrom_name} = $seq;
close(GENOME);
print "genome hashing completed\n";
#----------------------------------------------------------------

my $mRNA_counter = 0;
$line = <INPUT>; print OUT $line; print NO_MATCH $line;
while($line = <INPUT>)  { 
     chomp($line);
     @line = split /\t/,$line;
     if ($line[5] eq '++') {  print OUT "$line\n"; $mRNA_counter++;}
     else {            
            if (exists($annotation{$line[1]})) {
               $mRNA_counter++; 
               my $strand    = $annotation{$line[1]}[3];
               my $chrom_seq = $genome{$annotation{$line[1]}[0]};
                  if ($line[5] eq '-+')
                     { 
                        if($strand eq '+')
                         {  $new_5UTRstart = $annotation{$line[1]}[1] - 100 + $line[3];
                            $sequence5UTR = substr($chrom_seq, $new_5UTRstart - 1, $annotation{$line[1]}[1] - $new_5UTRstart);
                         }
                        else
                         {  $new_5UTRstart = $annotation{$line[1]}[2] + 100 - $line[3];
                            $sequence5UTR = substr($chrom_seq, $annotation{$line[1]}[2], $new_5UTRstart - $annotation{$line[1]}[2]);
                            $sequence5UTR = reverse($sequence5UTR);
                            $sequence5UTR =~ tr/ACTG/TGAC/;
                         }
                        print OUT "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$sequence5UTR"."$line[6]\n";
                     }
                  elsif ($line[5] eq '+-')
                     {
                        if($strand eq '+')
                         {  $new_3UTRstart = $annotation{$line[1]}[2] + 100 - $line[4];
                            $sequence3UTR = substr($chrom_seq, $annotation{$line[1]}[2], $new_3UTRstart - $annotation{$line[1]}[2]);
                         }
                        else
                         {  $new_3UTRstart = $annotation{$line[1]}[1] - 100 + $line[4];
                            $sequence3UTR = substr($chrom_seq, $new_3UTRstart - 1, $annotation{$line[1]}[1] - $new_3UTRstart);
                            $sequence3UTR = reverse($sequence3UTR);
                            $sequence3UTR =~ tr/ACTG/TGAC/;
                         }
                        print OUT "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$line[6]"."$sequence3UTR\n";
                     }
                  elsif ($line[5] eq '--')
                     {
                        if($strand eq '+')
                         {  $new_5UTRstart = $annotation{$line[1]}[1] - 100 + $line[3];
                            $new_3UTRstart = $annotation{$line[1]}[2] + 100 - $line[4];
                            $sequence5UTR = substr($chrom_seq, $new_5UTRstart - 1, $annotation{$line[1]}[1] - $new_5UTRstart);
                            $sequence3UTR = substr($chrom_seq, $annotation{$line[1]}[2], $new_3UTRstart - $annotation{$line[1]}[2]);
                         }
                        else
                         {  $new_5UTRstart = $annotation{$line[1]}[2] + 100 - $line[3];
                            $new_3UTRstart = $annotation{$line[1]}[1] - 100 + $line[4];
                            $sequence5UTR = substr($chrom_seq, $annotation{$line[1]}[2], $new_5UTRstart - $annotation{$line[1]}[2]);
                            $sequence3UTR = substr($chrom_seq, $new_3UTRstart - 1, $annotation{$line[1]}[1] - $new_3UTRstart);
                            $sequence5UTR = reverse($sequence5UTR);
                            $sequence3UTR = reverse($sequence3UTR);
                            $sequence5UTR =~ tr/ACTG/TGAC/;
                            $sequence3UTR =~ tr/ACTG/TGAC/;
                         }
                        print OUT "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$sequence5UTR"."$line[6]"."$sequence3UTR\n";
                     }
            }
            else { print NO_MATCH "$line\n";  }
     } 
  }
print "Number of mRNA transcripts:\t$mRNA_counter\n";
close (INPUT); close (OUT); close (NO_MATCH);


# Convert to 2-line fasta format
open (INPUT, "<mRNA_100"); open (OUT, ">mRNA_100.fasta");
$line = <INPUT>;
while ($line = <INPUT>) {
     @line = split /\t/, $line;
     print OUT ">$line[0]\n";
     print OUT $line[6];
}
close(INPUT); close(OUT);
unlink "mRNA_100"; unlink "no_Match"; unlink "temp3";















