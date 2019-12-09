# NextflowRNAseq

Pipeline NEXTFLOW for quality, trimming, indexation, STAR mapping, matrice generation with Cufflinks, RSEM and/or FEELnc.  

To run this nextflow pipeline, you have to :

1 / __Prepare your working directory__ : Decompress nextflow_star.tar.gz and copy your reads and annotation/genome files in specific directories.

2 / **Run two command lines to install nextflow** (module load bioinfo/Java8; curl -fsSL get.nextflow.io | bash) and to use the pipeline on Genologin (sbatch STARcuff.sh ...).


# A - PIPELINE DESCRIPTION:

## Choose a mode to choose a pipeline

This pipeline can launch several tools depending of the mode you choose to run. Several modes are availables for this Nextflow pipeline:

* QUANTIFICATION on a new model with FeatureCounts (--mode quantifFC): trimming - star - Cufflinks - Cuffmerge - FeatureCounts steps.
```console
cd /path/to/your/working/directory/; module purge; 
source /usr/local/bioinfo/src/MultiQC/venv_multiqc-v1.5/bin/activate; source /usr/local/bioinfo/src/cutadapt/cutadapt-1.14/venv_python-2.7.2/bin/activate; 
sbatch STARcuff.sh  --workdir /path/to/your/working/directory/ --resultsdir FC/ --genome Bos_taurus.UMD3.1.dna.toplevel.fa --annotation Bos_taurus.UMD3.1.92.gtf --conditions neg,pos --mode quantifFC; 
```

* QUANTIFICATION on a known model with FeatureCounts (--mode quantifFCRef): trimming - star - FeatureCounts steps.
```console
cd /path/to/your/working/directory/; module purge; 
source /usr/local/bioinfo/src/MultiQC/venv_multiqc-v1.5/bin/activate; source /usr/local/bioinfo/src/cutadapt/cutadapt-1.14/venv_python-2.7.2/bin/activate; 
sbatch STARcuff.sh  --workdir /path/to/your/working/directory/ --resultsdir FC/ --genome Bos_taurus.UMD3.1.dna.toplevel.fa --annotation Bos_taurus.UMD3.1.92.gtf --conditions neg,pos --mode quantifFCRef; 
```

* QUANTIFICATION on a known model with RSEM (--mode quantifRSEM) : trimming - star - RSEM steps 
```console
cd /path/to/your/working/directory/; module purge; 
source /usr/local/bioinfo/src/MultiQC/venv_multiqc-v1.5/bin/activate; source /usr/local/bioinfo/src/cutadapt/cutadapt-1.14/venv_python-2.7.2/bin/activate; 
sbatch STARcuff.sh  --workdir /path/to/your/working/directory/  --resultsdir TESTrsem/ --genome Bos_taurus.UMD3.1.dna.toplevel.fa --annotation Bos_taurus.UMD3.1.92.gtf --conditions neg,pos --mode quantifRSEM
```

* QUANTIFICATION with Cufflinks/Cuffmerge (--mode quantifCufflinks): trimming - star - Cufflinks - Cuffmerge steps
```console
cd /path/to/your/working/directory/; module purge; 
source /usr/local/bioinfo/src/MultiQC/venv_multiqc-v1.5/bin/activate; source /usr/local/bioinfo/src/cutadapt/cutadapt-1.14/venv_python-2.7.2/bin/activate; 
sbatch STARcuff.sh  --workdir /path/to/your/working/directory/ --resultsdir CUFF/ --genome Bos_taurus.UMD3.1.dna.toplevel.fa --annotation Bos_taurus.UMD3.1.92.gtf --conditions neg,pos --mode quantifCufflinks; 
```

* MODEL GENERATION with FEELnc and messagers classification (--mode model): trimming - star - cufflinks - cuffmerge - FEELnc steps
```console
cd /path/to/your/working/directory/; module purge; 
source /usr/local/bioinfo/src/MultiQC/venv_multiqc-v1.5/bin/activate; source /usr/local/bioinfo/src/cutadapt/cutadapt-1.14/venv_python-2.7.2/bin/activate; 
sbatch STARcuff.sh  --workdir /path/to/your/working/directory/  --resultsdir MODEL/ --genome Bos_taurus.UMD3.1.dna.toplevel.fa --annotation Bos_taurus.UMD3.1.92.gtf --conditions neg,pos --mode model
```

Warning: If you have no protein conding in your GTF reference or not enought protein coding (min 100), please add known_mRNA.gtf file in annotation/ directory, this file will be used instead of GTF reference.
 
## Description of each tool

* __Trimming adaptateurs__ with cutadapt version 1.8.3: 
```console
--minimum-length= 20
-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC 
-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
```

* __STAR Genome Generate__ (only once) 

* __STAR mapping__:
STAR run with paired read, cufflinks option, quantMode on transcriptome

Mapping parameters description <br/>
(source: http://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STAR.posix/doc/STARmanual.pdf):
```console
--readFilesCommand zcat
--readFilesCommand : Compressed inputs files
--outSAMtype : Output BAM sorted by coordonate
--outFileNamePrefix : Output files name prefix (including full or relative path)
--outFilterType BySJout : Keep only those reads that contain junctions that passed ltering into SJ.out.tab
--quantMode TranscriptomeSAM : Output SAM/BAM alignments to transcriptome into a separate file
--outSAMstrandField intronMotif : Strand derived from the intron motif. Reads with inconsistent and/or non-canonical introns are filtered out. ( Cuffinks-like strand field flag)
--outFilterIntronMotifs RemoveNoncanonical : Filter out alignments that contain non-canonical junctions.
--alignIntronMin 10 : Minimum intron size: genomic gap is considered intron if its length >=alignIntronMin, otherwise it is considered Deletion. 
--alignIntronMax 25000 : Maximum intron length 
--outFilterMismatchNmax 10 : Maximum number of mismatches per pair, large number switches off this filter.
```

* __Samtools sort__ aligned BAM files on names.

* __RSEM genome generate__ (only once), preparation of transcripts BAM files from STAR to be used by RSEM, then RSEM.
Performs gene and isoform level quantification from RNA-Seq data. RSEM is a software package that quantifies gene and isoform abundances from single-end (SE) or paired-end (PE) RNA-Seq data. The software enables visualization of its output through probabilistically-weighted read alignments and read depth plots. It does not require a reference genome and can be used for quantifying de novo transcriptome assemblies. 

* __Cufflinks__ v2.1.1 on merged BAM:
Cufflinks is both the name of a suite of tools and a program within that suite. Cufflinks the program assembles transcriptomes from RNA-Seq data and quantifies their expression. <br/>
(source: http://cole-trapnell-lab.github.io/cufflinks/cufflinks/index.html)

Parameters description :
```console
-g/GTF-guide reference_annotation.(gtf/gff) : Tells Cufflinks to use the supplied reference annotation a GFF file to guide RABT assembly. Reference transcripts will be tiled with faux-reads to provide additional information in assembly. Output will include all reference transcripts as well as any novel genes and isoforms that are assembled.
-F/min-isoform-fraction : After calculating isoform abundance for a gene, Cufflinks filters out transcripts that it believes are very low abundance, because isoforms expressed at extremely low levels often cannot reliably be assembled, and may even be artifacts of incompletely spliced precursors of processed transcripts. This parameter is also used to filter out introns that have far fewer spliced alignments supporting them. The default is 0.1, or 10% of the most abundant isoform (the major isoform) of the gene.
-j/pre-mrna-fraction : Some RNA-Seq protocols produce a significant amount of reads that originate from incompletely spliced transcripts, and these reads can confound the assembly of fully spliced mRNAs. Cufflinks uses this parameter to filter out alignments that lie within the intronic intervals implied by the spliced alignments. The minimum depth of coverage in the intronic region covered by the alignment is divided by the number of spliced reads, and if the result is lower than this parameter value, the intronic alignments are ignored. The default is 15%.
```

* __Cuffmerge__ transcripts from Cufflinks in a merged.gtf file.
Cuffmerge merge together several Cufflinks assemblies with a reference file.
When you have multiple RNA-Seq libraries and you’ve assembled transcriptomes from each of them, we recommend that you merge these assemblies into a master transcriptome. This step is required for a differential expression analysis of the new transcripts you’ve assembled. Cuffmerge performs this merge step.

* __FEELnc__ : 3 modules in one FEELnc script

1. FEELnc_filter.pl : __Extract, filter candidate transcripts__.

2. FEELnc_codpot.pl : __Compute the coding potential of candidate transcripts__.

The first step of the pipeline (FEELnc_filter) consists in filtering out unwanted/spurious transcripts and/or transcripts overlapping (in sense) exons of the reference annotation and especially protein_coding exons as they more probably correspond to new mRNA isoforms (see -b,--biotype option).<br/>
The main step of the pipeline (FEELnc_codpot) aims at computing the CPS i.e the coding potential score (between [0-1]) foreach of the candidate transcripts in the candidate_lncRNA.gtf file.<br/>

In the absence of species-specific lncRNAs set, machine-learning strategies require to simulate non-coding RNA sequences to train the model. We developed 2 approaches:

Mode __Shuffle__ : First approach involves that lncRNAs derived from "debris" of protein-coding genes (Duret et al. 2006). For this strategy that we called shuffle, the set of mRNAs are taken and shuffled while preserving 7-mer frequencies using Ushuffle. If you want to use the shuffle mode, please check that the fasta_ushuffle binary is in your PATH.

Mode __Intergenic__ : Another more naive approach called intergenic consists in extracting random sequences from the genome of interest to model species-specific noncoding sequences. In this case, the reference genome file is required (ref_genome.FA) and the mode of the lncRNA sequences simulation have to been set to intergenic.
As in the previous module, if your reference annotation file ("ref_annotation.GTF") contains additionnal fields such transcript_biotype and/or transcript_status in the GENCODE annotation or ENSEMBL, you can extract them manually or by using the -b option (as to get the best training set of known mRNAs.

3. FEELnc_classifier.pl: __Classify lncRNAs__ based on their genomic localization wrt others transcripts.

lncRNA.gtf = output file of step 2 = concatenat noORF gtf file and lncRNA gtf file because noORF are, by definition, lncRNA too.


* __MultiQC-v1.7__ (https://multiqc.info/):<br/>
Aggregate results from bioinformatics analyses across many samples into a single report
MultiQC searches a given directory for analysis logs and compiles a HTML report. It's a general use tool, perfect for summarising the output from numerous bioinformatics tools.

* Possible to add a __debug option__ "--debug Index,Cleaning" (debug steps to REMOVE: GenomeDir, Cleaning, STARmap, sort, prep, Rsem, Cuff, Mrg, FC, fref, FEELnc, Quality) if you want to run only one or several steps. If you don't need debug mode, just not add this parameter.<br/>
```console
GenomeDir      : no STAR Genome Dir
Cleaning       : no cutadapt
STARmap        : no STAR mapping
sort           : no samtools sort
prep           : no RSEM prepare reference
Rsem           : no RSEM
Cuff           : no Cufflinks
Mrg            : no Cuffmerge
FC             : no FeatureCounts on merged.gft generated by Cuffmerge
fref           : no FeatureCounts on reference GTF
FEELnc         : no 3 steps of FEELnc
Quality        : no multiQC and Nextflow statistics
```

# B - HOW TO USE NEXTFLOW SCRIPT:

## In order to run this pipeline in your computer you will required

    Unix-like operating system
    Java 7 (or higher)

## Your working directory in your working path :

Create 3 directories __annotation__ , __genome__ and __reads__ :

```console
  $ ls nextflow_star/
  annotation  genome  nextflow  nextflow.conf  reads  STARcuff_V2.sh  star_V2.nf README script_concat_rsem.py report_multiqc.sh
```

"genome" directory contains your fasta reference file.
"annotation" directory contains your GTF and GFF3 files.
"reads" directory contains all your samples files.

Add your sequences files in "reads" directory, your GTF file in "annotation" directory and your fasta reference file in "genome" directory.

```console
  $ ls genome/
  Bos_taurus.UMD3.1.dna.toplevel.fa

  $ ls annotation/
  Bos_taurus.UMD3.1.92.gtf  
                     
  $ ls reads/
  1-1-neg_CAGATC_L007_R1.fastq.gz
  1-1-neg_CAGATC_L007_R2.fastq.gz
  2-1-pos_CAGATC_L008_R1.fastq.gz
  2-1-pos_CAGATC_L008_R2.fastq.gz
  ...
```

__IMPORTANT__ : Reads files names *condition*_R1.fastq.gz and *condition*_R2.fastq.gz


Please check that all .sh files  are executables (-rwxr-xr-x) otherwise ...
```console
$ chmod a+x *.sh
```


## (Only if you don't run jobs on Genologin) - Install Nextflow entering the following command in the shell terminal


Nextflow is already installed on Genologin.

```console
$ module load bioinfo/Java8; curl -fsSL get.nextflow.io | bash
```

Then dependencies are downloaded. and the executable file have been created:
                                                               
      N E X T F L O W
      version 0.32.0 build 4897
      last modified 27-09-2018 10:17 UTC (12:17 CEST)
      cite doi:10.1038/nbt.3820
      http://nextflow.io

Nextflow installation completed. Please note:
- the executable file `nextflow` has been created in the folder: /home/Sarah/Documents/NEXTFLOW/TESTS/nextflow_star
- you may complete the installation by moving it to a directory in your $PATH

If you need more informations : User guide [instruction](https://www.nextflow.io/docs/latest/index.html)

## Before run pipeline, we have to specify some parameters


### 1 - In your work directory:


Please, copy your genome fasta file in /pathToYourWorkdir/genome/ directory and copy your GTF annotation file in /pathToYourWorkdir/annotation/ directory.
Moreover, check that your reads file are in fastq.gz format and paired with _R1 and _R2 names (see params.reads variable).

### 2 - Launch your pipeline execution using this command


First, connect to Genologin :

```console
  $ ssh -YX yourLogin@genologin.toulouse.inra.fr
```

Then, run Nextflow script (sbatch is specific to SLURM) :

```console
  $cd /path/to/YourWorkDir/; 
  $ source /usr/local/bioinfo/src/MultiQC/venv_multiqc-v1.5/bin/activate; source /usr/local/bioinfo/src/cutadapt/cutadapt-1.14/venv_python-2.7.2/bin/activate; module purge; 
  $ sbatch STARcuff.sh  --workdir /path/to/YourWorkDir/ --resultsdir resultsDirName --genome animal.dna.toplevel.fa --annotation animal.gtf --conditions condition1,condition2 --mode quantifFC
```

__Warnings__

* Please choose your mode quantifFC, quantifFCRef, quantifRSEM, quantifCufflinks, or model.

* "condition1" : if your sequences file name is, for instance, "reads/1-1-neg_CAGATC_L007_R1.fastq.gz", the first condition must be written "neg" in the following way (the expression "neg" must be found in the file name) . This Nextflow script is designed, for the time being, for two maximum conditions.

* Even if you run the same pipeline for multiple batches of samples, reference genome indexation will only be run once.

* slurm-JobNumber.out contains pipeline logs. If any problem, please share this file.

* A README file is generated in your results/ directory (/path/to/YourWorkDir/results/README). This file contains more informations about each step of the pipeline (inputs and outputs files and runned steps).

__Biblio__
https://www.biorxiv.org/content/10.1101/833962v1

