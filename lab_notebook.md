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

## Merging CDS fasta and repeat fasta files

```bash
cat urchin-families.fa GCF_000002235.5_Spur_5.0_rna.fna > combined_fasta.fa
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

## Kallisto running for bulk RNAseq quantification

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

## TEcandidates

```bash
TEcandidates.sh -t=32 -r=32 -c=1 -l=20 -te=Spur_repeats.gff3 -g=Spur_5.0.fa -fq=. -m=PE -N=10000 1> TEcandidates_urchin.log 2> TEcandidates_urchin.err
```

## Kallisto bustools 

```bash
bedtools getfasta -fi Spur_5.0.fa -bed expressed_RE.gff3 -split -s -name -fo expressed_repeats.fa
cat GCF_000002235.5_Spur_5.0_rna.fna expressed_repeats.fa > Genes+repeats.fa
kallisto index -i Genes_repeats Genes+repeats.fa

for i in SRR11599813 SRR11599814 SRR11599815 SRR11599816 SRR11599817 SRR11599818 SRR11599819 SRR11599820 SRR11599821 SRR11599822 SRR11599823 SRR11599824
do
  mkdir "$i"_out
  kb count -i Genes_repeats.idx -g tr2g.tsv -x 10xv3 -t 16 -o "$i"_out --mm $i* &> $i.log
done
```

