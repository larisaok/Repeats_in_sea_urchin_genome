# Laboratory notebook
## Downloading reference genome and annotation

```bash
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/235/GCF_000002235.5_Spur_5.0/GCF_000002235.5_Spur_5.0_genomic.fna.gz
gunzip GCF_000002235.5_Spur_5.0_genomic.fna
mv GCF_000002235.5_Spur_5.0_genomic.fna Spur.fna # renaming reference into Spur.fna
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/235/GCF_000002235.5_Spur_5.0/GCF_000002235.5_Spur_5.0_genomic.gff.gz
gunzip GCF_000002235.5_Spur_5.0_genomic.gff.gz

```

## RepeatModeler and RepeatMasker running

```bash
nohup BuildDatabase -name urchin Spur.fna >& building_database.out & # fast 
nohup RepeatModeler -pa 8 -database urchin >& run.out & # 96 hours on 8 CPU, 30 Gb memory
nohup RepeatMasker -pa 8 -lib urchin-families.fa -nolow -no_is Spur.fna & # 13 hours on 8 CPU, 30 Gb memory

```

## Getting gff file of repetitive elements

```bash
cd /opt/RepeatMasker
wget http://www.repeatmasker.org/utils/rmOutToGFF3.pl
chmod +x rmOutToGFF3.pl
cd ~/urchin
perl /opt/RepeatMasker/rmOutToGFF3.pl Spur.fna.out > transposons.gff
```

## Gff grooming

```bash
awk '{FS="\t"; OFS="\t"} {print}' transposons.gff | perl -pe 's/\ .*//' > clear_transposons.gff
awk '{FS="\t"; OFS="\t"} {print $1, $2, $9, $4, $5, $6, $7, $8, $3}' clear_transposons.gff > replased_transposons.gff
```

## Converting replased_transposon.gff to fasta

```bash
nohup bedtools getfasta -fi Spur.fna -bed replased_transposons.gff -name -s &> transposons_fasta &
```

## Cutting useless features from genome gff

```bash

sed '/tRNA/d' GCF_000002235.5_Spur_5.0_genomic.gff > cut_annot # get rid of tRNA
sed '/rRNA/d' cut_annot > cut_annot2 # get rid of rRNA
sed '/tandem_repeat/d' cut_annot2 > cut_annot3 # get rid of tandem repeats
awk -F"\t" '!seen[$4, $5]++' cut_annot3 > cut_annot4 # get rid of features with the same coordinates
```

## Converting gff to CDS multifasta

```bash
nohup bedtools getfasta -fi Spur.fna.masked -bed cut_annotation -name -s &> genes_fasta &
```

## Merging CDS fasta and repeat fasta files

```bash
cat transposons_fasta genes_fasta > combined_fasta.fa
```

## Downloading RNA-seq data (dataset1)

```bash
nano downloading_RNA-seq_data.sh
for i in SRR8863030 SRR8863031 SRR8863032 SRR8863033 SRR8863034 SRR8863035 SRR8863036 SRR8863037 SRR8863038 SRR8863039 SRR8863027 SRR8863028 SRR8863020
do
fastq-dump --split-3 $i
done
chmod +x downloading_RNA-seq_data.sh
nohup ./downloading_RNA-seq_data.sh &> downloading_RNA-seq_data.log &
```

## Quality control of RNA-seq data

```bash
mkdir fastqc_results
nano fastqc.sh
for i in SRR*
do
fastqc -o fastqc_results $i
done
chmod +x fastqc.sh
nohup ./fastqc.sh &> fastqc.log &
```

## Kallisto running

```bash
# index building (8 CPU, 120 Gb memory, took about 4 hours to complete)
nohup ../kallisto/kallisto index --make-unique -i all_transcripts.idx combined_fasta.fa &> kallisto_building_index_out &

for i in SRR88630*_1.fastq
do
  BASE=${i%%_1.fastq}
  mkdir ../kallisto_output/output_$BASE 
  kallisto quant -t 32 -i ../all_transcripts.idx -o ../kallisto_output/output_$BASE -b 30 $BASE"_1.fastq" $BASE"_2.fastq" &>> kallisto_running.log
  wait $(pgrep kallisto)
  echo
  echo
done
```
