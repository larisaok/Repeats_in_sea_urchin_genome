# Analysis of repetitive elements transcription during embryogenesis of sea urchin *Strongylocentrotus purpuratus*

## Introduction
Repetitive elements (RE) occupy significant part of eukaryotic genomes and some of them are shown to play diverse roles in genome regulation. For example, they are implied to be the origins of non-coding RNA. During embryogenesis of the sea urchin, a significant number of RE are expressed (up to 80%). The role of these elements in the regulation of biological processes remains unknown.

## Aim
The aim of this study was to identify the RE in the genome of the purple sea urchin and analyze their expression at different stages of embryogenesis.

## Project steps and methods
1. Identifing and classification RE in the genome *S. purpuratus* (Spur 5.0) using RepeatModeler
2. Quantify RE abundance in genome.
3. Reference transcriptome assembling which includes RE and previously annotated genes.
4. Quantification of transcripts abundances using Kallisto.
5. Statistical analysis in Sleuth.

Detailed description of steps 1-4 can be found in [lab notebook](https://github.com/larisaok/Repeats_in_sea_urchin_genome/blob/master/lab_notebook.md)
Step 5 description can be found in [sleuth_script](https://github.com/larisaok/Repeats_in_sea_urchin_genome/blob/master/sleuth_script.R)

## RNA-seq data

RNA-seq data of *S. purpuratus* embryos on different stages were obtained from SRA database (Project [PRJNA531297](https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP191285&o=acc_s%3Aa))

## System requirements
* RepeatModeler and RepeatMasker can be found [here](http://www.repeatmasker.org/RepeatModeler/)
* R version 3.5.3
* kallisto version 0.46.0 can be found [here](https://pachterlab.github.io/kallisto/starting)
* sleuth R package v0.30.0 ([conda package](https://anaconda.org/bioconda/r-sleuth))

## Results

For example,
![PCA](https://github.com/larisaok/Repeats_in_sea_urchin_genome/blob/master/pca.png) 
PCA plot illustrates that experimental points from embryos before and after gastrulation are grouped together, suggesting that pregastrulation stages shows higher similarity.

![Volcano plot](https://github.com/larisaok/Repeats_in_sea_urchin_genome/blob/master/volcano_plot.png) 
Volcano plot shows that majority of transcripts change their expression after gastrulation.

