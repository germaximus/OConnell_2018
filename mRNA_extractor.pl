#!/usr/bin/perl
# extracts all CDS from a given *.gbk file. Leaves 100 nt to the left and to the right.

###extract RefSeq records with gene names and CDS length from *.gbk file.
open (INPUT, "<$ARGV[0]"); open (OUT, ">temp1");
print OUT "gene_name\tRefSeq:ID\tCDS_length\t5UTR_length\t3UTR_length\tSequence\n";
while($line = <INPUT>)
   {
      @line=split /\s+/,$line;
      if (index($line[1], "NM_") != -1)
             {
                $accession = $line[1];
                until (($line[1] eq 'gene') && (scalar(@line) == 3) && (index($line[2], '..') != -1)) {$line = <INPUT>; @line=split /\s+/,$line;}
                chomp($line[2]);                                                                                       
                @gene_coordinate = split /\.\./,$line[2];
                until (index($line[1],'/gene="') != -1) {$line = <INPUT>; @line=split /\s+/,$line;}
                $gene = substr($line[1],7, -1);
                until (($line[1] eq 'CDS') && (scalar(@line) == 3)) {$line = <INPUT>; @line=split /\s+/,$line;}
                if(index($line[2],"join") == -1)  
                    {chomp($line[2]);
                     @CDS_coordinate = split /\.\./,$line[2];
                     $CDS_length = $CDS_coordinate[1] - $CDS_coordinate[0] + 1;
                     $UTR5_length = $CDS_coordinate[0] - 1;
                     $UTR3_length = $gene_coordinate[1] - $CDS_coordinate[1];
                     
                     until ($line[0] eq 'ORIGIN') {$line = <INPUT>; chomp ($line); @line=split /\s+/,$line;}
                     until (index($line, '//') != -1) {$line = <INPUT>;
                                                       chomp ($line);
                                                       @line=split /\s+/,$line;
                                                       shift(@line); shift(@line);
                                                       $sequence_line = join('',@line);
                                                       $sequence_line = uc($sequence_line);
                                                       $sequence_whole .= $sequence_line;}
                     print OUT "$gene\t$accession\t$CDS_length\t$UTR5_length\t$UTR3_length\t$sequence_whole\n"; }
                else {until (index($line, '//') != -1) {$line = <INPUT>;}}
                undef $sequence_whole;
             }
      else { until (index($line, '//') != -1) {$line = <INPUT>;} }
   }
close (INPUT); close (OUT);

###select the longest CDS's for every gene. If CDS are equal, the one with the longest 5'-UTR is chosen. If 5`-UTR are also equal, then 3'-UTR is evaluated.
open (INPUT, "<temp1"); open (OUT, ">temp2");
$line = <INPUT>; print OUT $line;
while($line = <INPUT>)
   {  chomp($line);
      @line = split ("\t", $line);
      @{$gene_name{$line[0]}{$line[1]}} = ($line[2], $line[3], $line[4], $line[5]);
   }
foreach $name (keys %gene_name)
   {
      my @CDS_lengths;
      my @UTR3_lengths;
      my @UTR5_lengths;
      my @CDSs;
      my @RefSeqID;
      my @array;
      my $CDS_counter=0;
      my $UTR5_counter=0;
      my $UTR3_counter=0;
      
      if (scalar keys%{$gene_name{$name}} == 1)
            {
               @RefSeqID = keys %{$gene_name{$name}};
               $accession = $RefSeqID[0];
               $CDS = $gene_name{$name}{$RefSeqID[0]}[0];
               $UTR5 = $gene_name{$name}{$RefSeqID[0]}[1];
               $UTR3 = $gene_name{$name}{$RefSeqID[0]}[2];
               $Sequence = $gene_name{$name}{$RefSeqID[0]}[3];
               print OUT "$name\t$accession\t$CDS\t$UTR5\t$UTR3\t$Sequence\n";
            }
      #some genes have more than one RefSeq ID, i.e. isoforms with varying 5' 3' UTRs or CDS lengths
      else {
               foreach $RefSeqID (keys %{$gene_name{$name}})
                    {  push(@CDS_lengths,$gene_name{$name}{$RefSeqID}[0]); } 
               @CDS_lengths = sort { $a <=> $b } @CDS_lengths; $CDS_longest = $CDS_lengths[-1];
               foreach $record (@CDS_lengths) {if($record eq $CDS_longest) {$CDS_counter++;}}
               if ($CDS_counter == 1)
                     {
                        foreach $RefSeqID (keys %{$gene_name{$name}})
                          { if($gene_name{$name}{$RefSeqID}[0] == $CDS_longest)
                             {
                               $accession = $RefSeqID;
                               $CDS = $gene_name{$name}{$RefSeqID}[0];
                               $UTR5 = $gene_name{$name}{$RefSeqID}[1];
                               $UTR3 = $gene_name{$name}{$RefSeqID}[2];
                               $Sequence = $gene_name{$name}{$RefSeqID}[3];
                               print OUT "$name\t$accession\t$CDS\t$UTR5\t$UTR3\t$Sequence\n";
                             }
                          }
                     }
               else  {  foreach $RefSeqID (keys %{$gene_name{$name}})
                           {if($gene_name{$name}{$RefSeqID}[0] == $CDS_longest) 
                           {push(@UTR5_lengths,$gene_name{$name}{$RefSeqID}[1]);}}
                        @UTR5_lengths = sort { $a <=> $b } @UTR5_lengths; $UTR5_longest = $UTR5_lengths[-1];
                        foreach $record (@UTR5_lengths) {if($record eq $UTR5_longest) {$UTR5_counter++;}}
                        if ($UTR5_counter == 1)
                              {
                                 foreach $RefSeqID (keys %{$gene_name{$name}})
                                   { if(($gene_name{$name}{$RefSeqID}[1] == $UTR5_longest) && ($gene_name{$name}{$RefSeqID}[0] == $CDS_longest))
                                      {
                                         $accession = $RefSeqID;
                                         $CDS = $gene_name{$name}{$RefSeqID}[0];
                                         $UTR5 = $gene_name{$name}{$RefSeqID}[1];
                                         $UTR3 = $gene_name{$name}{$RefSeqID}[2];
                                         $Sequence = $gene_name{$name}{$RefSeqID}[3];
                                         print OUT "$name\t$accession\t$CDS\t$UTR5\t$UTR3\t$Sequence\n";
                                      }
                                   }
                              }
                        else  {  foreach $RefSeqID (keys %{$gene_name{$name}})
                                   {if(($gene_name{$name}{$RefSeqID}[1] == $UTR5_longest) && ($gene_name{$name}{$RefSeqID}[0] == $CDS_longest))
                                   {push(@UTR3_lengths,$gene_name{$name}{$RefSeqID}[2]);}}
                                 @UTR3_lengths = sort { $a <=> $b } @UTR3_lengths; $UTR3_longest = $UTR3_lengths[-1];
                                 foreach $record (@UTR3_lengths) {if($record eq $UTR3_longest) {$UTR3_counter++;}}
                                   if ($UTR3_counter == 1)
                                        {
                                          foreach $RefSeqID (keys %{$gene_name{$name}})
                                            { if(($gene_name{$name}{$RefSeqID}[1] == $UTR5_longest) && ($gene_name{$name}{$RefSeqID}[0] == $CDS_longest) && ($gene_name{$name}{$RefSeqID}[2] == $UTR3_longest))
                                               {
                                                 $accession = $RefSeqID;
                                                 $CDS = $gene_name{$name}{$RefSeqID}[0];
                                                 $UTR5 = $gene_name{$name}{$RefSeqID}[1];
                                                 $UTR3 = $gene_name{$name}{$RefSeqID}[2];
                                                 $Sequence = $gene_name{$name}{$RefSeqID}[3];
                                                 print OUT "$name\t$accession\t$CDS\t$UTR5\t$UTR3\t$Sequence\n";
                                               }
                                            }
                                        }
                                    #some genes (for instance Tpm1) have isoforms with identical lengths of CDS, 5UTR and 3UTR. For them, I pick one record at random.    
                                    else { foreach $RefSeqID (keys %{$gene_name{$name}})
                                            { if(($gene_name{$name}{$RefSeqID}[1] == $UTR5_longest) && ($gene_name{$name}{$RefSeqID}[0] == $CDS_longest) && ($gene_name{$name}{$RefSeqID}[2] == $UTR3_longest))
                                               {
                                                 push(@array,$RefSeqID);
                                               }
                                            }
                                           $random_key = $array[0]; #not trully "random", just the first elemment of the array
                                           $accession = $random_key;
                                           $CDS = $gene_name{$name}{$random_key}[0];
                                           $UTR5 = $gene_name{$name}{$random_key}[1];
                                           $UTR3 = $gene_name{$name}{$random_key}[2];
                                           $Sequence = $gene_name{$name}{$random_key}[3];
                                           print OUT "$name\t$accession\t$CDS\t$UTR5\t$UTR3\t$Sequence\n";
                                         }
                              }
                     }
           }
   }
close (INPUT); close (OUT);

##extract CDS + 100 nt left/right from gene sequence. If 5UTR and/or 3UTR are shorter than 100 nt, raise a flag.
open (INPUT, "<temp2"); open (OUT, ">temp3");
$line = <INPUT>; @line=split /\t/,$line;
print OUT "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t"."Flag"."\t$line[5]";
while($line = <INPUT>)
  {      chomp ($line); @line=split /\t/,$line;
         if(($line[3] >= 100) && ($line[4] >= 100))
             {  $sequence = substr($line[5],$line[3] - 100, $line[2] + 200); $flag = '++';            }
         if(($line[3] >= 100) && ($line[4] < 100))
             {  $sequence = substr($line[5],$line[3] - 100, $line[2] + 100 + $line[4]); $flag = '+-'; }
         if(($line[3] < 100) && ($line[4] >= 100))
             {  $sequence = substr($line[5],0, $line[2] + $line[3] + 100); $flag = '-+';              }
         if(($line[3] < 100) && ($line[4] < 100))
             {  $sequence = substr($line[5],0, $line[2] + $line[3] + $line[4]); $flag = '--';         }
         print OUT "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$flag\t$sequence\n";   
  }
close (INPUT); close (OUT);
unlink (temp1, temp2);