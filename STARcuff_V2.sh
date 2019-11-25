#!/bin/bash
#SBATCH -p workq
#SBATCH -J nextflowSTAR
#SBATCH --time=1-01:00:00 #1 jour et 1h
#SBATCH  --mem 40G



args=("$@");
echo $# arguments passed;

# Le script doit recevoir au minimum 6 arguments
if [ $# -lt 6 ]
 then
          # Si le nombre d'arguments est inférieur à 6
          # on retourne le code erreur 1
          echo "Incorrect number of arguments";
          exit 128 #code erreur : invalide arguments number
  else
          # Si le nombre d'arguments est supérieur ou égal à 6
          # on retourne le code erreur 0
          echo "Correct number of arguments";
fi


# scaner programming env: these switches turn some bugs into errors
set -o errexit -o pipefail -o noclobber -o nounset

OPTIONS=mwgacr
LONGOPTS=mode,debug,workdir,resultsdir,genome,annotation,conditions

modewf=-
debugwf=PROD
wd=-
fasta=-
gtf=- 
cond=-
resultsdir=-

# now enjoy the options in order and nicely split until we see --
while [ $# != 0 ]; do
    case "$1" in
        -m|--mode)                  #modes : quantifRSEM, quantifFC, model, etc...
            modewf="$2"
            shift 2 
            ;;
        -d|--debug)                  #debug : GenomeDir, Cleaning, STARmap, prep, Rsem, Cuff, Mrg, FC, fref, FEELnc, Quality
            debugwf="$2"
            shift 2 
            ;;
        -w|--workdir)
            wd="$2"
            shift 2
            ;; 
        -r|--resultsdir)
            rd="$2"
            shift 2
            ;; 
        -g|--genome)
            fasta="$2"
            shift 2
            ;;
        -a|--annotation)
            gtf="$2"
            shift 2
            ;;
        -c|--conditions)
         #   shift
         #   cond="$@"
            cond="$2"
            shift 2
            ;;
        *)
            echo "bad parameters" 1>&2
            exit 1
            ;;
    esac
done

echo "workdir=$wd \n resultsdir=$rd \n mode=$modewf \n genome=$fasta \n annotation=$gtf \n debug=$debugwf \n";
#change neg,pos in neg pos
listcond=${cond//,/ }
echo "list of conditions = $listcond  \n";
#idem for debug
#debuglist=${debugwf//,/ }
#echo "Debug mode = $debuglist \n";

#smaman@genologin1 /work/smaman/nextflow_star_rsem_cufflinks_cuffmerge_feelnc $ ./STARcuff.sh  --workdir /work/smaman/nextflow_star_rsem_cufflinks_cuffmerge_feelnc/ --genome Bos_taurus.UMD3.1.dna.toplevel.fa --annotation Bos_taurus.UMD3.1.92.gtf --conditions neg,pos --mode quantif
#workdir=/work/smaman/nextflow_star_rsem_cufflinks_cuffmerge_feelnc/, mode=quantif, genome=Bos_taurus.UMD3.1.dna.toplevel.fa, annotation=-, conditions=neg,pos
#list of conditions = neg pos


#$wd= ${args[0]}
#$modewf
#$fasta=${args[1]}
#$gtf=${args[2]}
#$listcond=${args[3]}


#############################################  TESTS FILES #################################################################################################


######################## INPUTS FILES ######################################################################################################################
#TODO : Attention reads files names *condition*_R1.fastq.gz and *condition*_R2.fastq.gz

#smaman@genologin2 /work/smaman/nextflow_star_rsem_cufflinks_cuffmerge_feelnc/reads $ ls
#MT_rep1_1_Ch6.fastq.gz  MT_rep1_2_Ch6.fastq.gz  WT_rep1_1_Ch6.fastq.gz  WT_rep1_2_Ch6.fastq.gz
#smaman@genologin2 /work/smaman/nextflow_star_rsem_cufflinks_cuffmerge_feelnc/reads $ mv MT_rep1_1_Ch6.fastq.gz MT_rep1_Ch6_R1.fastq.gz     #NO _1 and _2
#smaman@genologin2 /work/smaman/nextflow_star_rsem_cufflinks_cuffmerge_feelnc/reads $ mv MT_rep1_2_Ch6.fastq.gz MT_rep1_Ch6_R2.fastq.gz
#smaman@genologin2 /work/smaman/nextflow_star_rsem_cufflinks_cuffmerge_feelnc/reads $ mv WT_rep1_1_Ch6.fastq.gz WT_rep1_Ch6_R1.fastq.gz
#smaman@genologin2 /work/smaman/nextflow_star_rsem_cufflinks_cuffmerge_feelnc/reads $ mv WT_rep1_2_Ch6.fastq.gz WT_rep1_Ch6_R2.fastq.gz

######################## ENV  #############################################################################################################################
#cd /work/smaman/nextflow_star_rsem_cufflinks_cuffmerge_feelnc/;  
#deactivate (si necessaire)
#source /usr/local/bioinfo/src/MultiQC/venv_multiqc-v1.5/bin/activate; source /usr/local/bioinfo/src/cutadapt/cutadapt-1.14/venv_python-2.7.2/bin/activate; module purge;


######################## 5 MODES       ####################################################################################################################

#FEATURECOUNTS + diminution du nombre de cpu (de 16 à 8)
#sbatch STARcuff.sh  --workdir /work/smaman/nextflow_star_rsem_cufflinks_cuffmerge_feelnc/ --resultsdir fc_tomate/ --genome ITAG2.3_genomic_Ch6.fasta --annotation ITAG_pre2.3_gene_models_Ch6.gtf --conditions MT,WT --mode quantifFC;
# ===> Job OK 
#tests with sh functions:
#sbatch STARcuff.sh  --workdir /work/smaman/nextflow_star_rsem_cufflinks_cuffmerge_feelnc/ --resultsdir fc_with_functions/ --genome ITAG2.3_genomic_Ch6.fasta --annotation ITAG_pre2.3_gene_models_Ch6.gtf --conditions MT,WT --mode quantifFC;
#sbatch STARcuff.sh  --workdir /work/smaman/nextflow_star_rsem_cufflinks_cuffmerge_feelnc/  --resultsdir fc_with_functions/ --genome ITAG2.3_genomic_Ch6.fasta --annotation ITAG_pre2.3_gene_models_Ch6.gtf --conditions MT,WT --mode quantifFC;


#FEATURECOUNTS REF   -  Diminution du nombre de cpu (de 16 à 8)
#sbatch STARcuff.sh  --workdir /work/smaman/nextflow_star_rsem_cufflinks_cuffmerge_feelnc/ --resultsdir fcREF_tomate/ --genome ITAG2.3_genomic_Ch6.fasta --annotation ITAG_pre2.3_gene_models_Ch6.gtf --conditions MT,WT --mode quantifFCRef;
#Nok / ref
# || Load annotation file ITAG_pre2.3_gene_models_Ch6.gtf ...                   ||
#  Failed to open the annotation file ITAG_pre2.3_gene_models_Ch6.gtf, or its format is incorrect, or it contains no 'exon' features.
# Expected format 
# SL2.40ch06	Cufflinks	exon	3688	4407	.	+	.	gene_id "XLOC_000001"; transcript_id "TCONS_00000001"; exon_number "1"; gene_name "Solyc06g005000.2.1"; oId "Solyc06g005000.2.1"; nearest_ref "Solyc06g005000.
#2.1"; class_code "="; tss_id "TSS1";
# Current format
#$ more annotation/ITAG_pre2.3_gene_models_Ch6.gtf
#SL2.40ch06	ITAG_eugene	exon	3688	4407	.	+	.	gene_id "Solyc06g005000.2.1"; transcript_id "Solyc06g005000.2.1";
#$cp CUFF_tomate/CUFFMERGE/merged.gtf annotation/.
#mv merged.gtf ITGA_merged.gtf
#sbatch STARcuff.sh  --workdir /work/smaman/nextflow_star_rsem_cufflinks_cuffmerge_feelnc/ --resultsdir fcREF_with_functions/ --genome ITAG2.3_genomic_Ch6.fasta --annotation ITGA_merged.gtf --conditions MT,WT --mode quantifFCRef;
#==> Job OK

#RSEM   - 8093150  log cutadapt OK avec multiqc OK + matrice de comptage mise en forme + Add -s "expected_count", "TPM", "FPKM"
#sbatch STARcuff.sh  --workdir /work/smaman/nextflow_star_rsem_cufflinks_cuffmerge_feelnc/ --resultsdir rsem_tomate/ --genome ITAG2.3_genomic_Ch6.fasta --annotation ITAG_pre2.3_gene_models_Ch6.gtf --conditions MT,WT --mode quantifRSEM;
# ===> Job OK
#tests with sh main and function: 
#sbatch STARcuff.sh  --workdir /work/smaman/nextflow_star_rsem_cufflinks_cuffmerge_feelnc/ --resultsdir rsem_with_functions/ --genome ITAG2.3_genomic_Ch6.fasta --annotation ITAG_pre2.3_gene_models_Ch6.gtf --conditions MT,WT --mode quantifRSEM;



#MODEL  FEELnc -   8085521  - with or without  known_mRNA.gtf in annotation directory  ( If you have no protein conding in your GTF reference or not enought protein coding (min 100), please add known_mRNA.gtf file in annotation/ directory, this file will be used instead of GTF reference).
#And add the possibility to have no {INPUT}.noORF.gtf file available
#sbatch STARcuff.sh  --workdir /work/smaman/nextflow_star_rsem_cufflinks_cuffmerge_feelnc/ --resultsdir MODEL_tomate/ --genome ITAG2.3_genomic_Ch6.fasta --annotation ITAG_pre2.3_gene_models_Ch6.gtf --conditions MT,WT --mode model;
#No protein_coding in ref GTF ....
#smaman@genologin2 /work/smaman/nextflow_star_rsem_cufflinks_cuffmerge_feelnc $  grep "protein_coding" /work/smaman/nextflow_star_rsem_cufflinks_cuffmerge_feelnc//annotation/ITAG_pre2.3_gene_models_Ch6.gtf
#smaman@genologin2 /work/smaman/nextflow_star_rsem_cufflinks_cuffmerge_feelnc $ 
#Add 101 differents protein_coding in ITGA gtf & 101 transcript_id different --> cp ITAG_pre2.3_gene_models_Ch6.gtf ITAG_pre2.3_gene_models_Ch6_PC.gtf
#Your input GTF file 'training_prot.gtf' contains only *1* transcripts !
#   Not enough to train the program (minimum == 100)...
# more ITAG_pre2.3_gene_models_Ch6.gtf
#SL2.40ch06	ITAG_eugene	exon	3688	4407	.	+	.	gene_id "protein_coding"; transcript_id "protein_coding";
#SL2.40ch06	ITAG_eugene	exon	3688	4407	.	+	.	gene_id "Solyc06g005000.2.1"; transcript_id "Solyc06g005000.2.1";
#sbatch STARcuff.sh  --workdir /work/smaman/nextflow_star_rsem_cufflinks_cuffmerge_feelnc/ --resultsdir MODEL_tomate/ --genome ITAG2.3_genomic_Ch6.fasta --annotation ITAG_pre2.3_gene_models_Ch6_PC.gtf --conditions MT,WT --mode model;
# ===> Job OK
#tests with sh functions:
#sbatch STARcuff.sh  --workdir /work/smaman/nextflow_star_rsem_cufflinks_cuffmerge_feelnc/ --resultsdir MODEL_with_functions/ --genome ITAG2.3_genomic_Ch6.fasta --annotation ITAG_pre2.3_gene_models_Ch6_PC.gtf --conditions MT,WT --mode model;


#CUFFLINKS
#sbatch STARcuff.sh  --workdir /work/smaman/nextflow_star_rsem_cufflinks_cuffmerge_feelnc/ --resultsdir CUFF_tomate/ --genome ITAG2.3_genomic_Ch6.fasta --annotation ITAG_pre2.3_gene_models_Ch6.gtf --conditions MT,WT --mode quantifCufflinks;
# ===> Job OK
#tests with sh functions:
#sbatch STARcuff.sh  --workdir /work/smaman/nextflow_star_rsem_cufflinks_cuffmerge_feelnc/  --resultsdir CUFF_with_functions/ --genome ITAG2.3_genomic_Ch6.fasta --annotation ITAG_pre2.3_gene_models_Ch6.gtf --conditions MT,WT --mode quantifCufflinks;


#CHANGES
#OK/ https://www.nextflow.io/docs/latest/tracing.html : -with-report $wd/$rd/starindex.report -with-timeline $wd/$rd/starindex.timeline  
# --> add in quality nextflow, each step is taken into account
#warning with DAG  : all processes are listed even if some of them are nor run. ( -with-dag flowchart.png )
# sbatch STARcuff.sh  --workdir /work/smaman/nextflow_star_rsem_cufflinks_cuffmerge_feelnc/ --resultsdir rsem_DAG/ --genome ITAG2.3_genomic_Ch6.fasta --annotation ITAG_pre2.3_gene_models_Ch6.gtf --conditions MT,WT --mode quantifRSEM; 
#to display file : cd /work/smaman/nextflow_star_rsem_cufflinks_cuffmerge_feelnc/rsem_with_functions; firefox Quantification_RSEM.timeline &

#OK / un seul multiqc : example: firefox rsem_with_functions/multiQC/multiqc_report.html &

#OK/ cut star_V2.nf script in several process to improve memory reservation. Example: add an optionnel process star and rsem dir.

#OK/ mv ${params.workpath}/${params.resultsdir}/Log.out ${params.multiQC}/Log.out;  --> done by mapping step in each mode so mv run in quality (last) step.

#debug:
#OK/debug (one or several steps excluded):a verif sur noFEELnc car des if inside (8131020)
#sbatch STARcuff.sh  --workdir /work/smaman/nextflow_star_rsem_cufflinks_cuffmerge_feelnc/ --debug noFEELnc --resultsdir MODEL_with_debug_nofeelnc/ --genome ITAG2.3_genomic_Ch6.fasta --annotation ITAG_pre2.3_gene_models_Ch6_PC.gtf --conditions MT,WT --mode model;
#OK/ test without 1 step:
#sbatch STARcuff.sh  --workdir /work/smaman/nextflow_star_rsem_cufflinks_cuffmerge_feelnc/ --resultsdir rsem_with_debug/ --debug noIndex --genome ITAG2.3_genomic_Ch6.fasta --annotation ITAG_pre2.3_gene_models_Ch6.gtf --conditions MT,WT --mode quantifRSEM;

#remove -resume overise some process are cached and not launch.

#TODO
#tests sur data perche cedric  

#best dag ?

#conditions parameters chanel ? cufflinks by condition , fcounts and rsem too ?


############################################################################################################################################################


export PERL5LIB=$PERL5LIB:/usr/local/bioinfo/src/FEELnc/FEELnc-v.0.1.1/lib/;


if [ $modewf == "quantifFC" ]
then
		if [ ! -d "$wd/$rd/" ]; then
		   mkdir $wd/$rd/; chmod 777 $wd/$rd/;
		fi
		echo "Your results directory is created : $wd/$rd/  \n";
		if [ ! -d "$wd/$rd/multiQC" ]; then
		   cd $wd/$rd/; mkdir multiQC/; chmod 777 multiQC/;
		fi
		echo "Your quality report  $wd/$rd/multiQC/ ";
 
        cd $wd/ &&
        
		#Step : STAR index
		#module load bioinfo/STAR-2.5.1b && module load bioinfo/samtools-1.4 &&
		./nextflow run star_V2.nf --mode starindex   --debug $debugwf --workpath $wd --resultsdir $rd  --fasta $fasta --gtf $gtf  && 
		#module unload bioinfo/STAR-2.5.1b && module unload bioinfo/samtools-1.4 &&
		echo "STAR index ....................................... FINISH \n" &&
		#Step : cutadapt then STAR map
		#module load bioinfo/cutadapt-1.14-python-3.4.3 && 
		./nextflow run star_V2.nf --mode starmapDir  --debug $debugwf --workpath $wd --resultsdir $rd   &&
		module load bioinfo/STAR-2.5.1b && module load bioinfo/samtools-1.4 &&
		./nextflow run star_V2.nf --mode starmap --debug $debugwf --workpath $wd  --resultsdir $rd  &&
		#module unload bioinfo/cutadapt-1.14-python-3.4.3 && 
		module unload bioinfo/STAR-2.5.1b && module unload bioinfo/samtools-1.4 &&
		echo "Cutadapt then STAR map ............................ FINISH \n" &&
		#Step : sort aligned to transcriptome file
		module load bioinfo/samtools-1.9 &&
		./nextflow run star_V2.nf --mode samtoolssort --debug $debugwf --workpath $wd --resultsdir $rd   &&
		module unload bioinfo/samtools-1.9 &&
		echo "Sort aligned to transcriptome files .............. FINISH \n" &&
		#Step : Cufflinks
		module load bioinfo/cufflinks-2.2.1 &&
		./nextflow run star_V2.nf --mode cufflinks --debug $debugwf --workpath $wd --resultsdir $rd --gtf $gtf --conditions  $listcond  &&
		module unload bioinfo/cufflinks-2.2.1 &&
		echo "Cufflinks  ........................................ FINISH \n" &&
		#Step : Cuffmerge
		module load bioinfo/cufflinks-2.2.1 &&
		./nextflow run star_V2.nf --mode cuffmerge --debug $debugwf --workpath $wd --resultsdir $rd --fasta $fasta --conditions  $listcond   &&
		module unload bioinfo/cufflinks-2.2.1 &&
		echo "Cuffmerge  ........................................ FINISH \n" &&
		#Step : featureCounts
		module load bioinfo/subread-1.6.0 &&
		./nextflow run star_V2.nf --mode fCounts  --debug $debugwf --workpath $wd --resultsdir $rd  &&
		module unload bioinfo/subread-1.6.0 &&
		echo "featureCounts .................................... FINISH \n" &&
		./nextflow run star_V2.nf --mode fCountsQuality  --debug $debugwf --workpath $wd --resultsdir $rd   &&
		./nextflow run star_V2.nf --mode qual --debug $debugwf --workpath $wd --resultsdir $rd -with-report $wd/$rd/Quantification_FeatureCounts.report -with-timeline $wd/$rd/Quantification_FeatureCounts.timeline  &&
		echo "multiQC .......................................... FINISH \n" &&
				
		echo "........................ YOUR NEXTFLOW PIPELINE IS FINISH \n" &&
		chmod 744 $wd/$rd/;
fi


if [ $modewf == "quantifFCRef" ]
then
		if [ ! -d "$wd/$rd/" ]; then
		   mkdir $wd/$rd/; chmod 777 $wd/$rd/;
		fi
		echo "Your results directory is created : $wd/$rd/  \n";
		if [ ! -d "$wd/$rd/multiQC" ]; then
		   cd $wd/$rd/; mkdir multiQC/; chmod 777 multiQC/;
		fi
		echo "Your quality report  $wd/$rd/multiQC/ ";
 
        cd $wd/ &&
        
       
		#Step : STAR index
		#module load bioinfo/STAR-2.5.1b && module load bioinfo/samtools-1.4 &&
		./nextflow run star_V2.nf --mode starindex  --debug $debugwf  --workpath $wd --resultsdir $rd  --fasta $fasta --gtf $gtf  && 
		#module unload bioinfo/STAR-2.5.1b && module unload bioinfo/samtools-1.4 &&
		echo "STAR index ....................................... FINISH \n" &&
		#Step : cutadapt then STAR map
		#module load bioinfo/cutadapt-1.14-python-3.4.3 && 
		./nextflow run star_V2.nf --mode starmapDir  --debug $debugwf --workpath $wd --resultsdir $rd   &&
		module load bioinfo/STAR-2.5.1b && module load bioinfo/samtools-1.4 &&
		./nextflow run star_V2.nf --mode starmap --debug $debugwf --workpath $wd  --resultsdir $rd  &&
		#module unload bioinfo/cutadapt-1.14-python-3.4.3 && 
		module unload bioinfo/STAR-2.5.1b && module unload bioinfo/samtools-1.4 &&
		echo "Cutadapt then STAR map ............................ FINISH \n" &&	
		#Step : sort aligned to transcriptome file
		module load bioinfo/samtools-1.9 &&
		./nextflow run star_V2.nf --mode samtoolssort --debug $debugwf --workpath $wd --resultsdir $rd   &&
		module unload bioinfo/samtools-1.9 &&
		echo "Sort aligned to transcriptome files .............. FINISH \n" &&
		#Step : featureCounts
		module load bioinfo/subread-1.6.0 &&
		./nextflow run star_V2.nf --mode fCountsOnRef --debug $debugwf --workpath $wd --resultsdir $rd --gtf $gtf  &&
		module unload bioinfo/subread-1.6.0 &&
		./nextflow run star_V2.nf --mode fCountsQualityOnRef  --debug $debugwf --workpath $wd --resultsdir $rd   &&
		echo "featureCounts .................................... FINISH \n" &&
		./nextflow run star_V2.nf --mode qual --debug $debugwf --workpath $wd --resultsdir $rd -with-report $wd/$rd/Quantification_FeatureCounts_OnRef.report -with-timeline $wd/$rd/Quantification_FeatureCounts_OnRef.timeline  &&
		echo "multiQC .......................................... FINISH \n" &&
				
		echo "........................ YOUR NEXTFLOW PIPELINE IS FINISH \n" &&
		chmod 744 $wd/$rd/;
fi


if [ $modewf == "quantifRSEM" ]
then
		if [ ! -d "$wd/$rd/" ]; then
		   mkdir $wd/$rd/; chmod 777 $wd/$rd/;
		fi
		echo "Your results directory is created : $wd/$rd/  \n";
		if [ ! -d "$wd/$rd/multiQC" ]; then
		   cd $wd/$rd/; mkdir multiQC/; chmod 777 multiQC/;
		fi
		echo "Your quality report  $wd/$rd/multiQC/ ";

        cd $wd/ &&
        
        
		#Step : STAR index 
		#module load bioinfo/STAR-2.5.1b && module load bioinfo/samtools-1.4 &&
     	./nextflow run star_V2.nf --mode starindex --debug $debugwf   --workpath $wd --resultsdir $rd  --fasta $fasta --gtf $gtf  && 
		#module unload bioinfo/STAR-2.5.1b && module unload bioinfo/samtools-1.4 &&
		echo "STAR index and RSEM Genome ........................ FINISH \n" &&
		#Step : cutadapt then STAR map
		./nextflow run star_V2.nf --mode starmapDir  --debug $debugwf --workpath $wd --resultsdir $rd   &&
		module load bioinfo/STAR-2.5.1b && module load bioinfo/samtools-1.4 &&
		./nextflow run star_V2.nf --mode starmap --debug $debugwf --workpath $wd --resultsdir $rd   &&
		module unload bioinfo/STAR-2.5.1b && module unload bioinfo/samtools-1.4 && 
		echo "Cutadapt then STAR map ............................ FINISH \n" &&
		#Step : sort aligned to transcriptome files
		module load bioinfo/samtools-1.9 &&
    	./nextflow run star_V2.nf --mode samtoolssort --debug $debugwf --workpath $wd --resultsdir $rd   &&
		module unload bioinfo/samtools-1.9 &&
		echo "Sort aligned to transcriptome files .............. FINISH \n" &&
		#Step :  RSEM 
		./nextflow run star_V2.nf --mode RSEMdir  --debug $debugwf --workpath $wd --resultsdir $rd   &&
		module load bioinfo/RSEM-1.3.0 &&
		./nextflow run star_V2.nf --mode RSEM --debug $debugwf --workpath $wd --resultsdir $rd --fasta $fasta --gtf $gtf  &&
		module unload bioinfo/RSEM-1.3.0 &&
		echo "RSEM .............................................. FINISH \n" &&
		#quality
		./nextflow run star_V2.nf --mode qual --debug $debugwf --workpath $wd --resultsdir $rd -with-dag $wd/$rd/RSEM_flowchart.png  -with-report $wd/$rd/Quantification_RSEM.report -with-timeline $wd/$rd/Quantification_RSEM.timeline  &&
		echo "multiQC .......................................... FINISH \n" &&
		
		
		echo "........................ YOUR NEXTFLOW PIPELINE IS FINISH \n" &&
		chmod 744 $wd/$rd/;
fi

if [ $modewf == "quantifCufflinks" ]
then
		if [ ! -d "$wd/$rd/" ]; then
		   mkdir $wd/$rd/; chmod 777 $wd/$rd/;
		fi
		echo "Your results directory is created : $wd/$rd/  \n";
		if [ ! -d "$wd/$rd/multiQC" ]; then
		   cd $wd/$rd/; mkdir multiQC/; chmod 777 multiQC/;
		fi
		echo "Your quality report  $wd/$rd/multiQC/ ";

        cd $wd/ &&
        
       
		#Step : STAR index 
		#module load bioinfo/STAR-2.5.1b && module load bioinfo/samtools-1.4 &&
		./nextflow run star_V2.nf --mode starindex  --debug $debugwf  --workpath $wd --resultsdir $rd  --fasta $fasta --gtf $gtf  && 
		#module unload bioinfo/STAR-2.5.1b && module unload bioinfo/samtools-1.4 &&
		echo "STAR index and RSEM Genome ........................ FINISH \n" &&
		
		#Step : cutadapt then STAR map
		./nextflow run star_V2.nf --mode starmapDir  --debug $debugwf --workpath $wd --resultsdir $rd   &&
		module load bioinfo/STAR-2.5.1b && module load bioinfo/samtools-1.4 &&
		./nextflow run star_V2.nf --mode starmap --debug $debugwf --workpath $wd --resultsdir $rd   &&
		module unload bioinfo/STAR-2.5.1b && module unload bioinfo/samtools-1.4 && 
		echo "Cutadapt then STAR map ............................ FINISH \n" &&
		#Step : sort aligned to transcriptome files
		module load bioinfo/samtools-1.9 &&
		./nextflow run star_V2.nf --mode samtoolssort --debug $debugwf --workpath $wd --resultsdir $rd   &&
		module unload bioinfo/samtools-1.9 &&
		echo "Sort aligned to transcriptome files .............. FINISH \n" &&
		#Step : Cufflinks
		module load bioinfo/cufflinks-2.2.1 &&
		./nextflow run star_V2.nf --mode cufflinks --debug $debugwf --workpath $wd --resultsdir $rd --gtf $gtf --conditions  $listcond  &&
		module unload bioinfo/cufflinks-2.2.1 &&
		echo "Cufflinks  ........................................ FINISH \n" &&
		#Step : Cuffmerge
		module load bioinfo/cufflinks-2.2.1 &&
		./nextflow run star_V2.nf --mode cuffmerge --debug $debugwf --workpath $wd --resultsdir $rd --fasta $fasta --conditions  $listcond   &&
		module unload bioinfo/cufflinks-2.2.1 &&
		echo "Cuffmerge  ........................................ FINISH \n" &&
		#quality
		./nextflow run star_V2.nf --mode qual --debug $debugwf --workpath $wd --resultsdir $rd  -with-dag $wd/$rd/Cufflinks_flowchart.png -with-report $wd/$rd/Quantification_Cufflinks.report -with-timeline $wd/$rd/Quantification_Cufflinks.timeline  &&
		echo "multiQC .......................................... FINISH \n" &&
			
		
		echo "......................... YOUR NEXTFLOW PIPELINE IS FINISH \n" &&
		chmod 744 $wd/$rd/;
fi


if [ $modewf == "model" ]
then
		if [ ! -d "$wd/$rd/" ]; then
		   mkdir $wd/$rd/; chmod 777 $wd/$rd/;
		fi
		echo "Your results directory is created : $wd/$rd/  \n";
		if [ ! -d "$wd/$rd/multiQC" ]; then
		   cd $wd/$rd/; mkdir multiQC/; chmod 777 multiQC/;
		fi
		echo "Your quality report  $wd/$rd/multiQC/ ";

        cd $wd/ &&
        
        #Step : STAR index and RSEM Genome 
		#module load bioinfo/STAR-2.5.1b && module load bioinfo/samtools-1.4 &&
		./nextflow run star_V2.nf --mode starindex  --debug $debugwf  --workpath $wd --resultsdir $rd  --fasta $fasta --gtf $gtf && 
		#module unload bioinfo/STAR-2.5.1b && module unload bioinfo/samtools-1.4 &&
		echo "STAR index and RSEM Genome ........................ FINISH \n" &&
		#Step : cutadapt then STAR map
		./nextflow run star_V2.nf --mode starmapDir  --debug $debugwf --workpath $wd --resultsdir $rd  &&
		module load bioinfo/STAR-2.5.1b && module load bioinfo/samtools-1.4 &&
		./nextflow run star_V2.nf --mode starmap  --debug $debugwf    --workpath $wd --resultsdir $rd   &&
		module unload bioinfo/STAR-2.5.1b && module unload bioinfo/samtools-1.4 && 
		echo "Cutadapt then STAR map ............................ FINISH \n" &&
		#Step : sort aligned to transcriptome files
		module load bioinfo/samtools-1.9 &&
		./nextflow run star_V2.nf --mode samtoolssort --debug $debugwf --workpath $wd --resultsdir $rd  &&
		module unload bioinfo/samtools-1.9 &&
		echo "Sort aligned to transcriptome files ............... FINISH \n" &&
		#Step : Cufflinks
		module load bioinfo/cufflinks-2.2.1 &&
		./nextflow run star_V2.nf --mode cufflinks --debug $debugwf --workpath $wd --resultsdir $rd --gtf $gtf --conditions  $listcond  &&
		module unload bioinfo/cufflinks-2.2.1 &&
		echo "Cufflinks  ........................................ FINISH \n" &&
		#Step : Cuffmerge
		module load bioinfo/cufflinks-2.2.1 &&
		./nextflow run star_V2.nf --mode cuffmerge --debug $debugwf --workpath $wd --resultsdir $rd --fasta $fasta --conditions  $listcond  &&
		module unload bioinfo/cufflinks-2.2.1 &&
		echo "Cuffmerge  ........................................ FINISH \n" &&
		#Step : Feelnc
		module load system/R-3.4.3 && module load bioinfo/FEELnc-v.0.1.1 && 
		#export PERL5LIB=$PERL5LIB:/usr/local/bioinfo/src/FEELnc/FEELnc-v.0.1.1/lib/ &&
		
		file="$wd/annotation/known_mRNA.gtf";
		if [[ -f "$file" ]];
		then
			echo "known_mRNA.gtf found.";
			./nextflow run star_V2.nf --mode FEELnc --debug $debugwf --workpath $wd --resultsdir $rd --fasta $fasta --gtf $wd/annotation/known_mRNA.gtf --conditions  $listcond ;
		else
			echo "known_mRNA.gtf not found, use $gtf.";
			./nextflow run star_V2.nf --mode FEELnc --debug $debugwf --workpath $wd --resultsdir $rd --fasta $fasta --gtf $gtf --conditions  $listcond ;
		fi

		module unload system/R-3.4.3 && module unload bioinfo/FEELnc-v.0.1.1 &&
		echo "Feelnc  ........................................... FINISH \n" &&
		./nextflow run star_V2.nf --mode qual --debug $debugwf --workpath $wd --resultsdir $rd  -with-dag $wd/$rd/MODEL_flowchart.png -with-report $wd/$rd/Model_Generation_with_FEELnc.report -with-timeline $wd/$rd/Model_Generation_with_FEELnc.timeline  &&
		echo "multiQC .......................................... FINISH \n" &&
		
		echo "......................... YOUR NEXTFLOW PIPELINE IS FINISH \n" &&
		chmod 744 $wd/$rd/;
fi
