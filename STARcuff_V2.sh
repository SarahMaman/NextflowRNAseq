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

#$wd= ${args[0]}
#$modewf
#$fasta=${args[1]}
#$gtf=${args[2]}
#$listcond=${args[3]}


export PERL5LIB=$PERL5LIB:/usr/local/bioinfo/src/FEELnc/FEELnc-v.0.1.1/lib/;


if [ $modewf == "quantifFC" ]
then
		if [ ! -d "$wd/$rd/" ]; then
		   mkdir $wd/$rd/; chmod 777 $wd/$rd/;
		fi
		echo "Your results directory is created : $wd/$rd/  \n";
		 
                cd $wd/ &&
                #Make output directories
                ./nextflow run star_V2.nf --mode quantifFCDir  --debug $debugwf --workpath $wd --resultsdir $rd   &&
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
		./nextflow run star_V2.nf --mode cufflinksListOK --debug $debugwf --workpath $wd --resultsdir $rd   &&
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
		#Step fatsQCreport
		module load bioinfo/FastQC_v0.11.5 &&
		./nextflow run star_V2.nf --mode QC  --debug $debugwf --workpath $wd --resultsdir $rd  &&
		echo "fastQCreport .................................... FINISH \n" &&		
				
		echo "........................ YOUR NEXTFLOW PIPELINE IS FINISH \n" &&
		chmod 744 $wd/$rd/;
fi


if [ $modewf == "quantifFCRef" ]
then
		if [ ! -d "$wd/$rd/" ]; then
		   mkdir $wd/$rd/; chmod 777 $wd/$rd/;
		fi
		echo "Your results directory is created : $wd/$rd/  \n";
		
                cd $wd/ &&
               #Make output directories
               ./nextflow run star_V2.nf --mode quantifFCRefDir  --debug $debugwf --workpath $wd --resultsdir $rd   &&
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
		#Step fatsQCreport
		module load bioinfo/FastQC_v0.11.5 &&
		./nextflow run star_V2.nf --mode QC  --debug $debugwf --workpath $wd --resultsdir $rd  &&
		echo "fastQCreport .................................... FINISH \n" &&
				
		echo "........................ YOUR NEXTFLOW PIPELINE IS FINISH \n" &&
		chmod 744 $wd/$rd/;
fi


if [ $modewf == "quantifRSEM" ]
then
		if [ ! -d "$wd/$rd/" ]; then
		   mkdir $wd/$rd/; chmod 777 $wd/$rd/;
		fi
		echo "Your results directory is created : $wd/$rd/  \n";
		cd $wd/ &&
                #Make output directiries 
                ./nextflow run star_V2.nf --mode quantifRSEMDir  --debug $debugwf --workpath $wd --resultsdir $rd   &&        
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
		./nextflow run star_V2.nf --mode qual --debug $debugwf --workpath $wd --resultsdir $rd  -with-report $wd/$rd/Quantification_RSEM.report -with-timeline $wd/$rd/Quantification_RSEM.timeline  &&
		echo "multiQC .......................................... FINISH \n" &&
		#Step fatsQCreport
		module load bioinfo/FastQC_v0.11.5 &&
		./nextflow run star_V2.nf --mode QC  --debug $debugwf --workpath $wd --resultsdir $rd  &&
		echo "fastQCreport .................................... FINISH \n" &&		
		
		echo "........................ YOUR NEXTFLOW PIPELINE IS FINISH \n" &&
		chmod 744 $wd/$rd/;
fi

if [ $modewf == "quantifCufflinks" ]
then
		if [ ! -d "$wd/$rd/" ]; then
		   mkdir $wd/$rd/; chmod 777 $wd/$rd/;
		fi
		echo "Your results directory is created : $wd/$rd/  \n";
                cd $wd/ &&
                #Make output directiries 
                ./nextflow run star_V2.nf --mode quantifCufflinksDir  --debug $debugwf --workpath $wd --resultsdir $rd   &&      
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
		./nextflow run star_V2.nf --mode cufflinksListOK --debug $debugwf --workpath $wd --resultsdir $rd   &&
		./nextflow run star_V2.nf --mode cuffmerge --debug $debugwf --workpath $wd --resultsdir $rd --fasta $fasta --conditions  $listcond   &&
		module unload bioinfo/cufflinks-2.2.1 &&
		echo "Cuffmerge  ........................................ FINISH \n" &&
		#quality
		./nextflow run star_V2.nf --mode qual --debug $debugwf --workpath $wd --resultsdir $rd  -with-report $wd/$rd/Quantification_Cufflinks.report -with-timeline $wd/$rd/Quantification_Cufflinks.timeline  &&
		echo "multiQC .......................................... FINISH \n" &&
		#Step fatsQCreport
		module load bioinfo/FastQC_v0.11.5 &&
		./nextflow run star_V2.nf --mode QC  --debug $debugwf --workpath $wd --resultsdir $rd  &&
		echo "fastQCreport .................................... FINISH \n" &&			
		
		echo "......................... YOUR NEXTFLOW PIPELINE IS FINISH \n" &&
		chmod 744 $wd/$rd/;
fi


if [ $modewf == "model" ]
then
		if [ ! -d "$wd/$rd/" ]; then
		   mkdir $wd/$rd/; chmod 777 $wd/$rd/;
		fi
		echo "Your results directory is created : $wd/$rd/  \n";
		
                cd $wd/ &&
               #Make output directiries 
               ./nextflow run star_V2.nf --mode modelDir  --debug $debugwf --workpath $wd --resultsdir $rd   &&
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
		./nextflow run star_V2.nf --mode cufflinksListOK --debug $debugwf --workpath $wd --resultsdir $rd   &&
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
		./nextflow run star_V2.nf --mode qual --debug $debugwf --workpath $wd --resultsdir $rd  -with-report $wd/$rd/Model_Generation_with_FEELnc.report -with-timeline $wd/$rd/Model_Generation_with_FEELnc.timeline  &&
		echo "multiQC .......................................... FINISH \n" &&
		#Step fatsQCreport
		module load bioinfo/FastQC_v0.11.5 &&
		./nextflow run star_V2.nf --mode QC  --debug $debugwf --workpath $wd --resultsdir $rd  &&
		echo "fastQCreport .................................... FINISH \n" &&
				
		echo "......................... YOUR NEXTFLOW PIPELINE IS FINISH \n" &&
		chmod 744 $wd/$rd/;
fi
