#!/usr/bin/perl
# WARNING!!! In this version I disabled interval based normalization. (scroll down for $end_row_sum).
# Usage:
# Coverage_processor.pl ARGV[0] ARGV[1]
# Where ARGV[0] - numerical, length of the CDSs used to select the subset, recommended 2000 for mice
#       ARGV[1] - binary "start" or "stop", specifies the location of the ribosome density profile
#       ARGV[2] - file name
# extract genes with specified length (centered either at start or stop codon)
open (FILE, "<$ARGV[2]"); open (OUT, ">temp1"); 
while ($line = <FILE>)
  {
     chomp($line);
     @line = split /\t/, $line;  
     if (scalar(@line) > $ARGV[0])
           {
             if($ARGV[1] eq "start")
              {  for(my $i = 0; $i <= $ARGV[0]; $i++)  {  print OUT "$line[$i]\t"; }
                 print OUT "\n";
              }
             else {die ("no start-stop argument is specified")}
           }
  }  close(FILE);close(OUT);

# calculate sum coverage for individual genes (row sum) and normalized average across each nucleotide for all genes (column average)
open (FILE, "<temp1"); open (OUT1, ">temp2"); open (OUT2, ">density_plot");
my @column_norm_full; my @column_norm_end;
my $counter = 0;
while ($line = <FILE>)
   {
     chomp($line); $line =~ s/\s+$//;
     @line = split /\t/, $line;
     $full_row_sum = 0;  # gene normalization factor based on the average coverage across the entire ORF
     #$end_row_sum = 0; # gene normalization factor based on the average cosumverage across the [ARGV[0] - 400nt] interval
             
     for (my $i = 1; $i < scalar(@line); $i++)   { $full_row_sum += $line[$i]; #if (scalar(@line) - $i <= 400) {$end_row_sum += $line[$i];}
                                                 }
     $full_average = $full_row_sum / (scalar(@line)-1);
     #$end_average = $end_row_sum / 400; # I use 400 nt from the end of the specified gene length to calculate the normalization coefficient
   
     if($full_average >= 0.25) # value 0.25 is arbitrary and designed to cut off low covered gene
       {    $counter++;
            my(@row_elements_norm_full, @row_elements_norm_end);
            for (my $i = 1; $i < scalar(@line); $i++)
             {
               $row_elements_norm_full[$i] = $line[$i] / $full_average;
               #$row_elements_norm_end[$i] = $line[$i] / $end_average;
               $column_norm_full[$i] += $row_elements_norm_full[$i];
               #$column_norm_end[$i] += $row_elements_norm_end[$i];
             }
            print OUT1 "$line[0]"."_full\t";
            shift (@row_elements_norm_full);  foreach $element (@row_elements_norm_full) {print OUT1 "$element\t";}
            print OUT1 "$full_average\n";
            #print OUT1 "$line[0]"."_end\t";
            #shift (@row_elements_norm_end);  foreach $element (@row_elements_norm_end) {print OUT1 "$element\t";}
            #print OUT1 "$end_average\n";
       }
   }
shift (@column_norm_full); shift (@column_norm_end);
for (my $i = 0; $i < scalar(@column_norm_full); $i++)
  {  $density_full_norm = $column_norm_full[$i] / $counter;
    # $density_end_norm =  $column_norm_end[$i] / $counter;
     print OUT2 "$density_full_norm\n";
    # print OUT2 "$density_end_norm\n";
  }  undef $counter;
close(FILE);close(OUT1);close(OUT2);
unlink temp1; unlink temp2; # these files were only created for manual troubleshooting
