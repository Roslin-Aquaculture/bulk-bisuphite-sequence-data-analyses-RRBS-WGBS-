# bulk-bisuphite sequence data analyses (RRBS and WGBS) by Robert Mukiibi and Diego Robledo
This is a pipeline that can used for analysing reduced representation bisulphite sequence (RRBS) and whole genome bisulphite sequence data (WGBS).
Here we incoporates the **bismark tools pipile** which are available here [FelixKrueger](https://github.com/FelixKrueger/Bismark) or [here](https://github.com/FelixKrueger/Bismark/archive/0.22.3.tar.gz) for bisulphite sequence data processing, and edgeR statitical analysis to for differential methylation of the profiled methylation sites.
Guide to using bismak tools is available [here](https://rawgit.com/FelixKrueger/Bismark/master/Docs/Bismark_User_Guide.html).
Differential methylation analyses is followed up by annotation of the differentially methylated CpG sites for biological or functional interpretation.
This pipeline has been utilized in our recently published article [Mukiibi et al.(2022)](https://doi.org/10.1016/j.ygeno.2022.110503).
**Note:** We are assuming paired-end dataset in the pipeline.

## 1. Bisulphite sequence data quality assessment using FastQC
[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) provides some basic sequencing quality on the raw reads including length distribution of the reads,per base and sequence quality scores, adapter content etc. These are useful for quality control in data cleaning prior to downstream analyses or processing.
```
outputFolder = /path_to_output_folder/
fastqc sample_R1.fastq -o ${outputFolder}
```
###### 2. Cleaning sequence data of low quality data using Trim Galore 
[Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) trims off adapters and sequence bases that of low quality scores. It also removes very short reads that might arise after quality processing.
```
trim_galore --paired --fastqc sample_R1.fastq sample_R2.fastq
```
###### 3. Insilico reference genome conversion (genome preparation) using bismark_genome_preparation
bismark_genome_preparation coverts reference genome to C->T converted genome version and a G->A version.
bismark_The converted genome versions are them indexed by [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) or [HISAT2](http://daehwankimlab.github.io/hisat2/) 
```
path_to_bowtie2=/usr/bin/bowtie2/ #give the full path of where your bowtie2 excecutable file is located
path_to_reference_genome=/path_to_folder_containing_your_reference_genome/ #give a full path to the folder where your reference genome is located.
bismark_genome_preparation --path_to_aligner ${path_to_bowtie2} --verbose ${path_to_reference_genome}
```
###### 4. Read alignment using bismark
bismark aligns the sequence read data to the converted reference genomes
bismak works [samtools](https://sourceforge.net/projects/samtools/files/samtools/) this will help choosing the alignment output file format *.bam or *.sam, without loading samtools the output will be sam file which is far heavy
```
samtools_path=/path_to_samtools/ #give the path to where you samtools is located
bismark  --bam --genome /path_to_folder_containing_your_reference_genome/ --samtools_path -1 sample1_1_val_1.fq -2 sample1_2_val_2.fq
```
###### 5. Methylation extraction using bismark_methylation_extractor
bismark_methylation_extractor takes the alignment file from bismark and extracts methylation call for each C analysed from this file.
For each sample bismark_methylation_extractor outputs a coverage file (*.cov) that contains all the CpG profiled, their methylation calls and coverage. 
It is these *.cov files that are utilized in the proceeding step of differential methylation analysis.
```
bismark_methylation_extractor --bedGraph --gzip sample1_1_val_1_bismark_bt2_pe.bam.gz
```
