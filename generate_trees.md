# Codes for phylogenetic tree reconstruction of functional genes from shotgun sequenced metagenomes-

### Note: You can skip Step1 and Step3 by running Step2 directly on nt.sequences followed by Step4. Or Skip Step3 and run Step4 directly after Step2.

## Step1- Get protein sequences from nt.sequences-
```
module load emboss/6.6.0 #name of the software

mkdir proteinsequences
transeq -table 11 -sequence annotated-gene1.fasta -outseq proteinsequences/protein-annotated-gene1.fasta #transeq is the name of the command used.
```

## Step2- Get alignment of protein sequences-
```
module load mafft/7.305 #name of the software

cd proteinsequences
mkdir proteinalignment
mafft --maxiterate 1000 --localpair protein-annotated-gene1.fasta > proteinalignment/aligned-protein-annotated-gene1.fasta #mafft is the name of the command used.
```

## Step3- Get protein coding nucleotide sequences-
```
cd proteinalignment
mkdir protein_nucleotides
pal2nal.pl aligned-protein-annotated-gene1.fasta ../../annotated-gene1.fasta -codontable 11 -output fasta > protein_nucleotides/protein_nucleo-annotated-gene1.fasta

##NOTE- If there is an error at this step, check out (A) of EXTRA
#pal2nal.pl is the name of the software used.
```

## Step4- Run IQTREE to get phylogenetic trees-
```
iqtree -s protein_nucleo-annotated-gene1.fasta -nt AUTO -alrt 2000 -bb 2000 -rcluster 10 -m TESTNEW  -pre gene1 -st CODON11

#-nt AUTO : automatically detect the number of nodes to use
#SH-aLRT (alrt) and UFBoot (bb) supports. One would typically start to rely on the clade if its SH-aLRT >= 80% and UFboot >= 95%
#-s ->to specify the alignment file
#output files- iqtree file, treefile, log file
#st CODON11 ->for protein coding genes of bacteria. You can use GTR+I+G also for the preliminary analysis.
```


# EXTRA-

## (A) If the protein sequences contain stop codons (asterisk in between the protein sequence)
```
#if the sequences contain stop codons in the middle-
#remove the stop codons from protein sequence file-
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' protein-alltermites-GTDBdb-COG0099.fasta | sed 's/.$//' | awk -F '\t'  '!($2 ~ /\*/)' | tr "\t" "\n" > nostop-protein-alltermites-GTDBdb-COG0099.fasta

#get the same headers out for nt sequence file-
>for i in nostop-protein-alltermites-GTDBdb-COG0*;do grep ">" ${i} > header-${i}.txt;done
 [output]
  229-01--HSRNA1--Archotermopsidae--Hodotermopsis-sjosdedti--China-Paleo-arctic--JPLABEBE_01158_d__Bacteria_p__Bacteroidota_c__Bacteroidia_1
  229-01--HSRNA1--Archotermopsidae--Hodotermopsis-sjosdedti--China-Paleo-arctic-- JPLABEBE_02417_d__Bacteria_p__Elusimicrobiota_c__Endomicrobia_o__Endomicrobiales_f__Endomicrobiaceae_1

>for i in header-nostop-protein-alltermites-GTDBdb-COG0*;do cat ${i} | sed 's/>//g' | sed 's/..$//' > 2-${i};done
  [output]
  229-01--HSRNA1--Archotermopsidae--Hodotermopsis-sjosdedti--China-Paleo-arctic--JPLABEBE_01158_d__Bacteria_p__Bacteroidota_c__Bacteroidia
  229-01--HSRNA1--Archotermopsidae--Hodotermopsis-sjosdedti--China-Paleo-arctic--JPLABEBE_02417_d__Bacteria_p__Elusimicrobiota_c__Endomicrobia_o__Endomicrobiales_f__Endomicrobiaceae

#get the nt. sequences out-
>seqtk subseq alltermites-GTDBdb-COG0185.fasta proteinsequences/2-header-nostop-protein-alltermites-GTDBdb-COG0185.fasta.txt > nostop-alltermites-GTDBdb-COG0185.fasta

#if the no. of sequences in ntsequence != protein sequence, then convert the nostop nt sequence to protein and start again-
module load emboss/6.6.0
for i in nostop-alltermites-GTDBdb-COG0*;do transeq -table 11 -sequence ${i} -outseq proteinsequences/protein-${i};done

```

## (B) If you want to do length specific cutoff-
```
## get lengths of sequences-
module load emboss/6.6.0
for i in annotated*fasta;do infoseq -auto -only -name -length ${i} > info-${i}.txt;done

## get mean length-
for i in info-*;do awk '{x+=$2; next} END{print x/NR}' ${i} > meanlength-${i}; done
cat meanlength* > meanlengthCOGs.txt

## get 70% of mean length per COG-
ls . > filenames.txt
paste -d ',' filenames.txt meanlengthCOGs.txt > meanlengthTreponemaCOGs.txt
awk -F"," '{print $0","70*$2/100}' meanlengthTreponemaCOGs.txt > 2-meanlengthTreponemaCOGs.txt #70% of mean length in column3

## extract the sequences based on length cutoff-
while read line;do column1=`echo ${line}| awk -F"," '{print $1}'`; column3=`echo ${line}| awk -F"," '{print $3}'`; perl extract500bpsabove.pl ${column3} Treponema-${column1}.fasta > lengthcutoff/lengthcut-${column1}.fasta;done < 2-meanlengthTreponemaCOGs.txt

#"extract500bpsabove.pl" perl script is attached in the Github repository
````
