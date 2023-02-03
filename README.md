# bulk-bisuphite sequence data analyses (RRBS and WGBS) by Robert Mukiibi
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
## 2. Cleaning sequence data of low quality data using Trim Galore 
[Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) trims off adapters and sequence bases that of low quality scores. It also removes very short reads that might arise after quality processing.
```
trim_galore --paired --fastqc sample_R1.fastq sample_R2.fastq
```
## 3. Insilico reference genome conversion (genome preparation) using bismark_genome_preparation
bismark_genome_preparation coverts reference genome to C->T converted genome version and a G->A version.
bismark_The converted genome versions are them indexed by [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) or [HISAT2](http://daehwankimlab.github.io/hisat2/) 
```
path_to_bowtie2=/usr/bin/bowtie2/ #give the full path of where your bowtie2 excecutable file is located
path_to_reference_genome=/path_to_folder_containing_your_reference_genome/ #give a full path to the folder where your reference genome is located.
bismark_genome_preparation --path_to_aligner ${path_to_bowtie2} --verbose ${path_to_reference_genome}
```
## 4. Read alignment using bismark
bismark aligns the sequence read data to the converted reference genomes
bismak works [samtools](https://sourceforge.net/projects/samtools/files/samtools/) this will help choosing the alignment output file format *.bam or *.sam, without loading samtools the output will be sam file which is far heavy
```
samtools_path=/path_to_samtools/ #give the path to where you samtools is located
bismark  --bam --genome /path_to_folder_containing_your_reference_genome/ --samtools_path -1 sample1_1_val_1.fq -2 sample1_2_val_2.fq
```
## 5. Methylation extraction using bismark_methylation_extractor
bismark_methylation_extractor takes the alignment file from bismark and extracts methylation call for each C analysed from this file.
For each sample bismark_methylation_extractor outputs a coverage file (*.cov) that contains all the CpG profiled, their methylation calls and coverage. 
It is these *.cov files that are utilized in the proceeding step of differential methylation analysis.
```
bismark_methylation_extractor --bedGraph --gzip sample1_1_val_1_bismark_bt2_pe.bam.gz
```
## 6. Differential methylation analysis using the R package edgeR 
In this step we utilize the Bioconductor package edgeR to analyse the the coverage data following the manual by [Chen et al.2017](https://f1000research.com/articles/6-2055/v2).
```
targets = read.csv("sample_files_and_descriptionFile.csv", header = T)##this file contains the sample COVERAGE(.cov)file name and the treatment or sample description
files <- paste0(targets$File,"")
yall <- readBismark2DGE(files, sample.names=targets$Sample)
chromosomes = read.csv("genomic_chromosomes_genes.csv", header = T)##some genomes have mitochondrial genome + scafold that are not mapped to chromosomes/ usually these are removed by providing a vector of chromosomes via this file
yall$samples$group <- factor(targets$Group)
table(yall$genes$Chr)
yall<-yall[yall$genes$Chr%in%chromosomes$Scafold_chr,]
##total up methylated and unmethylated "C"s to get total coverage.
Methylation = gl(2,1,ncol(yall), labels = c("Me", "Un"))
Me = yall$counts[, Methylation=="Me"]
Un = yall$counts[, Methylation=="Un"]
Coverage = Me + Un
##Filter our CpGs with less than 8 counts per sample.
HasCoverage = rowSums(Coverage >= 8) == length(targets$Sample) ## reamove any CpG with less that 8 counts per  sample
HasBoth=rowSums(Me) > 0 & rowSums(Un) > 0 ##Keep CpGs that show variability across samples, not fixed methylation or unmethylation state.
table(HasCoverage,HasBoth)
y = yall[HasCoverage & HasBoth,, keep.lib.sizes=F]
TotalLibSize = y$samples$lib.size[Methylation=="Me"] + y$samples$lib.size[Methylation=="Un"]
y$samples$lib.size = rep(TotalLibSize, each = 2)
Me = y$counts[, Methylation == "Me"]
Un = y$counts[, Methylation == "Un"]
#set GLM model
designSL = model.matrix(~0+Group, targets)
design = modelMatrixMeth(designSL)
y <- estimateDisp(y, design, trend="none")
y$common.dispersion
fit <- glmFit(y, design)
#Run contrast
contr = makeContrasts(GroupTreatment1vsTreatment2 = GroupTreatment1 - GroupTreatment2, levels = design) ###make contrast between two groups under investigation.
lrt = glmLRT(fit, contrast = contr)
summary(decideTests(lrt))
lrt_top <- topTags(lrt,adjust.method = "BH", sort.by = "PValue")
FDR_THRESHOLD = 0.05
dgegenes<- lrt_top$table[ lrt_top$table$FDR< FDR_THRESHOLD, ]
write.csv(dgegenes, file = "differentially_Methylated_CpGs.csv")###Save the differentially methylated CpGs
```
