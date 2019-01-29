# Guidance for preparing promoter-enhancer interactions
```sh
wget ftp://ccg.vital-it.ch/epdnew/H_sapiens/current/Hs_EPDnew.bed
perl -ne 's/ /\t/g;print' Hs_EPDnew.bed|cut -f1-4 > Hs_EPDnew.hg38.bed
perl -ne 'chomp;@t=split(/\t/);print "$t[0]\t",$t[1]<500000?1:$t[1]-500000,"\t",$t[2]+500000,"\t$t[3]\n"' Hs_EPDnew.hg38.bed > Hs_EPDnew.hg38.1Mb.bed
#Batch Coordinate Conversion (liftOver) using https://genome.ucsc.edu/cgi-bin/hgLiftOver
perl -ne '@t=split(/\t/);print "$t[0]\t",$t[1]-1000,"\t",$t[2]+1000,"\t$t[3]"' ~/Downloads/hglft_genome_4df9_ce7190.bed > Hs_EPDnew.hg19.1kb.bed
perl -ne '@t=split(/\t/);print "$t[0]\t",$t[1]-1000,"\t",$t[2]+1000,"\t$t[3]"' Hs_EPDnew.hg38.bed > Hs_EPDnew.hg38.1kb.bed
perl -ne '@t=split(/\t/);print "$t[0]\t$t[3]\t$t[4]\t$1\n" if $t[8]=~/genehancer_id=(\w+)/' genehancer46.csv | sort -k1,1 -k2,2n > genehancer46.hg38.bed
#Batch Coordinate Conversion (liftOver) using https://genome.ucsc.edu/cgi-bin/hgLiftOver
mv ~/Downloads/hglft_genome_41ab_ceacc0.bed genehancer46.hg19.bed
cat Hs_EPDnew.hg19.1kb.bed genehancer46.hg19.bed |perl -ne '@t=split(/\t/);print if $t[2]-$t[1]>100' | sort -k1,1 -k2,2n > promenh.hg19.bed
cat Hs_EPDnew.hg38.1kb.bed genehancer46.hg38.bed |perl -ne '@t=split(/\t/);print if $t[2]-$t[1]>100' | sort -k1,1 -k2,2n > promenh.hg38.bed
perl -ne 'chomp;@t=split(/\t/);print "$t[0]\t",$t[1]<499000?1:$t[1]-499000,"\t",$t[2]+499000,"\t$t[3]\n"' Hs_EPDnew.hg19.1kb.bed > Hs_EPDnew.hg19.1Mb.bed
bedtools intersect -a Hs_EPDnew.hg19.1Mb.bed -b ../primary_cohort_TAD_boundaries/H1.IS.All_TAD.bed -wo -f 0.6 -loj|perl -ne 'chomp;@t=split(/\t/);if ($t[4] eq "."){print "$t[0]\t$t[1]\t$t[2]\t$t[3]\n";next}print "$t[0]\t",$t[1]>$t[5]?$t[1]:$t[5],"\t",$t[2]<$t[6]?$t[2]:$t[6],"\t$t[3]\n"' > Hs_EPDnew.hg19.H1.1Mb.bed
bedtools intersect -a promenh.hg19.bed -b Hs_EPDnew.hg19.H1.1Mb.bed -wo -f 0.5 |perl -ne '@t=split(/\t/);print "$t[3]\t$t[7]\t",sprintf("%d",($t[1]+$t[2]-$t[5]-$t[6])/2),"\n" unless $t[3]=~/_\d+$/' > H1.hg19.pair
perl -ne 'chomp;@t=split(/\t/);for($i=-250;$i<=250;$i++){print "$t[0]\t",$t[1]+2000*($i-0.5),"\t",$t[1]+2000*($i+0.5),"\t$t[3]\t$i\n" if $t[1]+2000*($i-0.5)>0}' ~/Downloads/hglft_genome_1457_b49730.bed > Hs_EPDnew.hg19.win.bed
perl -ne '$i++;chomp;print "$_\t$i\n"' hg19.win2k.bed > hg19.win2k.bedi
bedtools intersect -a Hs_EPDnew.hg19.win.bed -b hg19.win2k.bedi -f 0.51 -wo |cut -f4,5,9 > Hs_EPDnew.hg19.win.idx
perl -ne 'chomp;@t=split(/\t/);print "$s\n$_\t" if $t[0] ne $old;$s=$_;$old=$t[0]' Hs_EPDnew.hg19.win.idx|cut -f1-3,6 > Hs_EPDnew.hg19.win.idx2
perl jsonProm.pl Hs_EPDnew.hg19.win.idx2 > genes.json
bedtools intersect -a genehancer46.hg19.bed -b hg19.win2k.bedi -wo |sort -k9,9rn|perl -ne '@t=split(/\t/);print "$t[3]\t$t[7]\n" unless defined $h{$t[3]};$h{$t[3]}=1' > genehancer46.hg19.idx
perl -ne 'print "chr,id\n" unless defined $old;@t=split(/\t/);print "$t[0],$t[3]" if $t[0] ne $old;$old=$t[0]' hg19.win2k.bedi > chrId.csv
#create intra-TAD P-E distances
bedtools intersect -a Hs_EPDnew.hg19.1Mb.bed -b ../primary_cohort_TAD_boundaries/H1.IS.All_TAD.bed  -wo -f 0.6 -loj|perl -ne 'chomp;@t=split(/\t/);if ($t[4] eq "."){print "$t[0]\t$t[1]\t$t[2]\t$t[3]\n";next}print "$t[0]\t",$t[1]>$t[5]?$t[1]:$t[5],"\t",$t[2]<$t[6]?$t[2]:$t[6],"\t$t[3]\n"' > Hs_EPDnew.hg19.H1.1Mb.bed
bedtools intersect -a Hs_EPDnew.hg19.1Mb.bed -b ../primary_cohort_TAD_boundaries/MES.IS.All_TAD.bed -wo -f 0.6 -loj|perl -ne 'chomp;@t=split(/\t/);if ($t[4] eq "."){print "$t[0]\t$t[1]\t$t[2]\t$t[3]\n";next}print "$t[0]\t",$t[1]>$t[5]?$t[1]:$t[5],"\t",$t[2]<$t[6]?$t[2]:$t[6],"\t$t[3]\n"' > Hs_EPDnew.hg19.MES.1Mb.bed
bedtools intersect -a Hs_EPDnew.hg19.1Mb.bed -b ../primary_cohort_TAD_boundaries/MSC.IS.All_TAD.bed -wo -f 0.6 -loj|perl -ne 'chomp;@t=split(/\t/);if ($t[4] eq "."){print "$t[0]\t$t[1]\t$t[2]\t$t[3]\n";next}print "$t[0]\t",$t[1]>$t[5]?$t[1]:$t[5],"\t",$t[2]<$t[6]?$t[2]:$t[6],"\t$t[3]\n"' > Hs_EPDnew.hg19.MSC.1Mb.bed
bedtools intersect -a Hs_EPDnew.hg19.1Mb.bed -b ../primary_cohort_TAD_boundaries/NPC.IS.All_TAD.bed -wo -f 0.6 -loj|perl -ne 'chomp;@t=split(/\t/);if ($t[4] eq "."){print "$t[0]\t$t[1]\t$t[2]\t$t[3]\n";next}print "$t[0]\t",$t[1]>$t[5]?$t[1]:$t[5],"\t",$t[2]<$t[6]?$t[2]:$t[6],"\t$t[3]\n"' > Hs_EPDnew.hg19.NPC.1Mb.bed
bedtools intersect -a Hs_EPDnew.hg19.1Mb.bed -b ../primary_cohort_TAD_boundaries/TRO.IS.All_TAD.bed -wo -f 0.6 -loj|perl -ne 'chomp;@t=split(/\t/);if ($t[4] eq "."){print "$t[0]\t$t[1]\t$t[2]\t$t[3]\n";next}print "$t[0]\t",$t[1]>$t[5]?$t[1]:$t[5],"\t",$t[2]<$t[6]?$t[2]:$t[6],"\t$t[3]\n"' > Hs_EPDnew.hg19.TPC.1Mb.bed
bedtools intersect -a Hs_EPDnew.hg19.1Mb.bed -b ../primary_cohort_TAD_boundaries/HC.IS.All_TAD.bed  -wo -f 0.6 -loj|perl -ne 'chomp;@t=split(/\t/);if ($t[4] eq "."){print "$t[0]\t$t[1]\t$t[2]\t$t[3]\n";next}print "$t[0]\t",$t[1]>$t[5]?$t[1]:$t[5],"\t",$t[2]<$t[6]?$t[2]:$t[6],"\t$t[3]\n"' > Hs_EPDnew.hg19.HC.1Mb.bed
bedtools intersect -a promenh.hg19.bed -b Hs_EPDnew.hg19.H1.1Mb.bed  -wo -f 0.5 |perl -ne '@t=split(/\t/);print "$t[3]\t$t[7]\t",($t[6]-$t[1]-500000),"\n" unless $t[3]=~/_\d+$/' > H1.hg19.pair
bedtools intersect -a promenh.hg19.bed -b Hs_EPDnew.hg19.HC.1Mb.bed  -wo -f 0.5 |perl -ne '@t=split(/\t/);print "$t[3]\t$t[7]\t",($t[6]-$t[1]-500000),"\n" unless $t[3]=~/_\d+$/' > HC.hg19.pair
bedtools intersect -a promenh.hg19.bed -b Hs_EPDnew.hg19.MES.1Mb.bed -wo -f 0.5 |perl -ne '@t=split(/\t/);print "$t[3]\t$t[7]\t",($t[6]-$t[1]-500000),"\n" unless $t[3]=~/_\d+$/' > MES.hg19.pair
bedtools intersect -a promenh.hg19.bed -b Hs_EPDnew.hg19.MSC.1Mb.bed -wo -f 0.5 |perl -ne '@t=split(/\t/);print "$t[3]\t$t[7]\t",($t[6]-$t[1]-500000),"\n" unless $t[3]=~/_\d+$/' > MSC.hg19.pair
bedtools intersect -a promenh.hg19.bed -b Hs_EPDnew.hg19.NPC.1Mb.bed -wo -f 0.5 |perl -ne '@t=split(/\t/);print "$t[3]\t$t[7]\t",($t[6]-$t[1]-500000),"\n" unless $t[3]=~/_\d+$/' > NPC.hg19.pair
bedtools intersect -a promenh.hg19.bed -b Hs_EPDnew.hg19.TPC.1Mb.bed -wo -f 0.5 |perl -ne '@t=split(/\t/);print "$t[3]\t$t[7]\t",($t[6]-$t[1]-500000),"\n" unless $t[3]=~/_\d+$/' > TPC.hg19.pair
```