#!/usr/bin/perl
# Usage: BLASTNprocesor.pl blast_result.txt

open (FILE, "<$ARGV[0]"); open (OUT1, ">redundant.temp"); open (OUT2, ">uniq.temp");
my %hash;
my %check_hash;
while(my $read=<FILE>) {
     my @read = split /\t/,$read;
     my $name = $read[0];
     push(@{$hash{$name}}, $read);
}
foreach my $name (sort keys %hash)   {
     my $counter = 0;
     if(scalar(@{$hash{$name}}) == 1 ) {    print OUT2 $hash{$name}[0];
                                            $check_hash{$name} = "unique";    }
     
     else { for (my $i = 1; $i < scalar(@{$hash{$name}}); $i++)
               {
                    my @subelements = split /\t/, $hash{$name}[$i];
                    if(($subelements[2] > 90) && ($subelements[3] > 50))
                        {
                          $check_hash{$subelements[0]} = "multiple blast hits"; #caution, it rewrites check_hash value if it already exists
                          $check_hash{$subelements[1]} = "multiple blast hits";
                          $counter++;
                        }
                }
            if ($counter == 0)
                {
                    print OUT2 $hash{$name}[0];
                    $check_hash{$name} = "unique";
                }                
            else
                {
                    foreach $element (@{$hash{$name}}) {print OUT1 $element;}
                }
          }
}
close(FILE); close(OUT1); close(OUT2);

#write down the fna file with unique non-redundant gene sequences
my $counter = 0;
my %refrence;
open (REF, "<mRNA_100.fasta");
while ($name = <REF>) {
    chomp($name);
    $name = substr($name, 1);
    $seq = <REF>;
    chomp($seq);
    $reference{$name} = $seq;
}
print "Total number of transcripts:\t";
print scalar(keys %reference)."\n";
close(REF);


open (FILE, "<uniq.temp"); open (OUT, ">mRNA_100uniq.fasta");
while (my $line = <FILE>) {
    my @line = split /\t/,$line;
    if(exists($reference{$line[0]})) {print OUT ">$line[0]\n"; print OUT $reference{$line[0]}."\n"; $counter++;}    
}
print "Number of unique transcripts:\t$counter\n";
close(FILE); close(OUT);
unlink("uniq.temp", "redundant.temp");


















