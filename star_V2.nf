params.workpath        = ''
params.debug           = ''
params.resultsdir      = ''
params.conditions      = ''
params.name            ="RNAseq_pipelines"
params.fasta           =''
params.gtf             =''
params.genome          ="${params.workpath}/genome/${params.fasta}"
/*
*params.annotation      ="${params.workpath}/annotation/${params.gtf}"
*/
params.annotation      ="${params.workpath}/annotation/*.gtf"
params.annotationGFF3  ="${params.workpath}/annotation/*.gff3"
params.reads           ="${params.workpath}/reads/*_{R1,R2}.fastq.gz"
params.readsPath       ="${params.workpath}/reads/"
params.readstrimmedDir ="${params.workpath}/${params.resultsdir}/TRIMMING/"
params.output          ="${params.workpath}/${params.resultsdir}/"
params.fastQCreport    ="${params.workpath}/${params.resultsdir}/fastQC_REPORT"
params.mappingDir      ="${params.workpath}/${params.resultsdir}/STAR_mapping" 
params.genomeDir       ="${params.workpath}/${params.resultsdir}/STAR_GenomeDir" 
params.rsemDir         ="${params.workpath}/${params.resultsdir}/RSEM_GenomeDir" 
params.cufflinksDir    ="${params.workpath}/${params.resultsdir}/CUFFLINKS"  
params.cuffmergeDir    ="${params.workpath}/${params.resultsdir}/CUFFMERGE"
params.feelnc          ="${params.workpath}/${params.resultsdir}/FEELNC"
params.rsem            ="${params.workpath}/${params.resultsdir}/RSEM"
params.fCounts         ="${params.workpath}/${params.resultsdir}/FeatureCounts"
params.fCountsOnRef    ="${params.workpath}/${params.resultsdir}/FeatureCountsOnRef"
params.multiQC         ="${params.workpath}/${params.resultsdir}/multiQC"
params.adaptatora      ="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
params.adaptA          ="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"
params.mode            = '' 


/*
 * Sigenae _ Sarah Maman _ 05-10-2018
 */


log.info "Trimming STAR Cufflinks - N F  ~  version 0.1"
log.info "============================================="
log.info "name                   : ${params.name}"
log.info "debug                  : ${params.debug}"
log.info "genome                 : ${params.genome}"
log.info "reads                  : ${params.reads}"
log.info "reads path             : ${params.readsPath}"
log.info "resultsdir             : ${params.resultsdir}"
log.info "annotation             : ${params.annotation}"
log.info "output                 : ${params.output}"
log.info "mode                   : ${params.mode}" 
log.info "quality                : ${params.multiQC}" 
log.info "\n"


/*
 * Creates the `read_pairs` channel that emits for each read-pair a tuple containing
 * three elements: the pair ID, the first read-pair file and the second read-pair file
 */
Channel
    .fromFilePairs( params.reads )                                             
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }  
    .set { read_pairs } 


/*
 * All the files ending with the suffix .fstq.gz are sent over the channel names seq. 
 * Then, these files are received by the process which will execute a fastqc report query on each of them.
 * seq = Channel.fromPath('${params.workpath}/reads/*.fastq.gz')
 */
seq = Channel.fromPath('reads/*.fastq.gz')
    
process indexation_star {

    cpus 8
    memory '60 GB'
	executor 'SLURM'
	
	publishDir "${params.genomeDir}"

	input:
	file 'params.genome' 
	file 'params.annotation' 
	
	
    script:
    //
    //STAR indexation
    //
    if( params.mode == 'starindex' )
    """         
		if [[ "${params.debug}" =~ .*GenomeDir.*  ]]; then
			echo "\n OK - NO STAR Genome Generate" >> ${params.output}/README;      
		else        
			if [ ! -d "${params.genomeDir}" ]; then
				cd ${params.workpath}/${params.resultsdir}; mkdir STAR_GenomeDir; chmod 777  STAR_GenomeDir; 
			fi  
			cd ${params.workpath}/${params.resultsdir};
	        module load bioinfo/STAR-2.5.1b; module load bioinfo/samtools-1.4; 
			STAR --runThreadN 8 --runMode genomeGenerate --limitGenomeGenerateRAM 35000000000 --genomeDir STAR_GenomeDir/ --genomeFastaFiles ${params.genome} --sjdbGTFfile  ${params.annotation} --sjdbOverhang 100;
			module unload bioinfo/STAR-2.5.1b; module unload bioinfo/samtools-1.4;
			echo "\n OK - STAR Genome Generate - Run at  `date`" >> ${params.output}/README;
		fi
    """
    else 
    """
		echo "No STAR indexation. \n";
    """
    }

process star_directories {

    executor 'SLURM'
	
    script:
    //
    //Directories creation
    //
    if( params.mode == 'starmapDir' )
    """
		if [ ! -d "${params.readstrimmedDir}" ]; then
			cd ${params.workpath}/${params.resultsdir}; mkdir  ${params.readstrimmedDir}/; chmod 777 ${params.readstrimmedDir}/;
		fi	
		if [ ! -d "${params.mappingDir}" ]; then
			cd ${params.workpath}/${params.resultsdir}; mkdir  ${params.mappingDir}/; chmod 777 ${params.mappingDir}/;
		fi	
    """
    else 
    """
        echo "Directories creation. \n";
    """
    }



process RSEM_directories {

    executor 'SLURM'
	
    script:
    //
    //Directories creation
    //  
	if( params.mode == 'RSEMdir' )
    """ 
		if [ ! -d "${params.rsemDir}" ]; then
			cd ${params.workpath}/${params.resultsdir}; mkdir  ${params.rsemDir}; chmod 777  ${params.rsemDir};
		fi
		if [ ! -d "${params.rsem}" ]; then
			cd ${params.workpath}/${params.resultsdir}; mkdir ${params.rsem}/; chmod 777  ${params.rsem}/; 
		fi
    """
    else 
    """
		echo "Directories creation. \n";
    """
    }



process star {
    
    cpus 8
    memory '60 GB'
	executor 'SLURM'
	

	input:
	file 'params.genome' 
	file 'params.annotation'
	set pair_id, file(reads) from read_pairs
	

    script:
    //
    //star mapping, samtools sort, cufflinks
    //
    if( params.mode == '' )
    """ 
		echo "\n No mode have been chosen.`" >> ${params.output}/README;
    """   
	else if( params.mode == 'starmap' )
    """
		cd ${params.workpath}/${params.resultsdir};
		if [ ! -d "${params.readstrimmedDir}/${pair_id}/" ]; then
			cd  ${params.readstrimmedDir}/; mkdir ${pair_id}/; chmod 777 ${pair_id}/; 
		fi
		cd ${params.readsPath};
		if [[ "${params.debug}" =~ .*Cleaning.*  ]]; then
			echo "\n OK - NO Cleaning" >> ${params.output}/README;      
		else  
			cutadapt -a ${params.adaptatora} -A ${params.adaptA} --minimum-length 20 --output ${params.readstrimmedDir}/${pair_id}/trimmed1.fastq --paired-output ${params.readstrimmedDir}/${pair_id}/trimmed2.fastq  ${reads}  > ${params.readstrimmedDir}/${pair_id}/cutadapt.log;
			echo "\n OK - Trimming on ${reads} - Run at  `date`" >> ${params.output}/README;
		fi
		cp ${params.readstrimmedDir}/${pair_id}/cutadapt.log ${params.multiQC}/${pair_id}_cutadapt.log;
		
		if [ ! -d "${params.mappingDir}/${pair_id}/" ]; then
			cd  ${params.mappingDir}/; mkdir ${pair_id}/; chmod 777 ${pair_id}/; 
		fi	
       
		if [[ "${params.debug}" =~ .*STARmap.*  ]]; then
			echo "\n OK - NO STAR mapping" >> ${params.output}/README;      
		else 
			STAR --genomeDir ${params.genomeDir} --readFilesIn ${params.readstrimmedDir}/${pair_id}/trimmed1.fastq ${params.readstrimmedDir}/${pair_id}/trimmed2.fastq --runThreadN 8  --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${params.mappingDir}/${pair_id}/ -outFilterType BySJout --quantMode TranscriptomeSAM --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --alignIntronMin 10 --alignIntronMax 25000 --outFilterMismatchNmax 10 ;
			echo "\n OK - STAR mapping on ${params.readstrimmedDir}/${pair_id}/trimmed1.fastq   and ${params.readstrimmedDir}/${pair_id}/trimmed2.fastq  - Run at  `date`" >> ${params.output}/README;

			cp ${params.mappingDir}/${pair_id}/SJ.out.tab ${params.multiQC}/${pair_id}_SJ.out.tab;
			cp ${params.mappingDir}/${pair_id}/Log.progress.out ${params.multiQC}/${pair_id}_Log.progress.out;
			cp ${params.mappingDir}/${pair_id}/Log.final.out ${params.multiQC}/${pair_id}_Log.final.out;
			cp ${params.mappingDir}/${pair_id}/Log.out ${params.multiQC}/${pair_id}_Log.out;
		fi
	"""
	else if( params.mode == 'samtoolssort' )
	"""
		if [[ "${params.debug}" =~ .*sort.*  ]]; then
			echo "\n OK - NO samtools sort" >> ${params.output}/README;      
		else
			cd ${params.workpath}/${params.resultsdir};
			chmod -R 777  ${params.mappingDir}/*;
			
			if  [ ! -d "${params.mappingDir}/${pair_id}/Aligned.toTranscriptome.sort.bam" ]; then
				samtools sort -n -@8 -o  ${params.mappingDir}/${pair_id}/Aligned.toTranscriptome.sort.bam  ${params.mappingDir}/${pair_id}/Aligned.toTranscriptome.out.bam;
			fi
			echo "\n OK - Samtools sort on ${pair_id} - Run at  `date`" >> ${params.output}/README;	
		fi	
    """     
    else 
    """
		echo "STAR step not yet or finish. \n";
    """
    }

Channel
    .fromFilePairs( params.reads )                                             
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }  
    .set { read_pairs2 } 
    
process Rsem {
    
    cpus 8
	executor 'SLURM'

	input:
	file 'params.genome' 
	file 'params.annotation'
	set pair_id, file(reads) from read_pairs2

    script:
    if( params.mode == 'RSEM' )
    """
		cd ${params.workpath}/${params.resultsdir};
	    
		if [ ! -d "${params.rsem}/${pair_id}/" ]; then
			cd  ${params.rsem}/; mkdir ${pair_id}/; chmod 777  ${pair_id}/; 
		fi
        
        if [[ "${params.debug}" =~ .*prep.*  ]]; then
			echo "\n OK - NO RSEM reference preparation" >> ${params.output}/README;      
        else 
			#module load bioinfo/RSEM-1.3.0; 
			rsem-prepare-reference --gtf ${params.annotation}  ${params.genome} ${params.rsemDir}/animal_refseq;
			#module unload bioinfo/RSEM-1.3.0; 
			echo "\n OK - RSEM Genome Generate - Run at  `date`" >> ${params.output}/README;
		fi
		
		if [[ "${params.debug}" =~ .*Rsem.*  ]]; then
            echo "\n OK - NO RSEM" >> ${params.output}/README;      
        else 	
			#Prepare BAM for RSEM
			convert-sam-for-rsem  ${params.mappingDir}/${pair_id}/Aligned.toTranscriptome.sort.bam ${params.rsem}/${pair_id}/Aligned.toTranscriptome.sort.bam.ConvertedBam;
			#RSEM 
			rsem-calculate-expression --bam --no-bam-output --estimate-rspd --calc-ci --seed 12345 -p 4 --ci-memory 30000 --paired-end ${params.rsem}/${pair_id}/Aligned.toTranscriptome.sort.bam.ConvertedBam.bam  ${params.rsemDir}/animal_refseq ${params.rsem}/${pair_id}/Aligned.toTranscriptome.out.bam.ConvertedBam_Quant;
		   
			cp ${params.rsem}/${pair_id}/Aligned.toTranscriptome.out.bam.ConvertedBam_Quant.stat/Aligned.toTranscriptome.out.bam.ConvertedBam_Quant.cnt   ${params.multiQC}/${pair_id}_Aligned.toTranscriptome.out.bam.ConvertedBam_Quant.cnt;
			cp ${params.rsem}/${pair_id}/Aligned.toTranscriptome.out.bam.ConvertedBam_Quant.stat/Aligned.toTranscriptome.out.bam.ConvertedBam_Quant.model  ${params.multiQC}/${pair_id}_Aligned.toTranscriptome.out.bam.ConvertedBam_Quant.model;
			cp ${params.rsem}/${pair_id}/Aligned.toTranscriptome.out.bam.ConvertedBam_Quant.stat/Aligned.toTranscriptome.out.bam.ConvertedBam_Quant.theta  ${params.multiQC}/${pair_id}_Aligned.toTranscriptome.out.bam.ConvertedBam_Quant.theta;
			
			#Script de concatenation des fichiers de quantifications des isoforms en matrice globale 
			python ${params.workpath}/script_concat_rsem.py -o ${params.rsem}/MATRICE_isoforms_expected_count -s expected_count -i ${params.rsem}/*/*isoforms*;
			python ${params.workpath}/script_concat_rsem.py -o ${params.rsem}/MATRICE_isoforms_TPM -s TPM -i ${params.rsem}/*/*isoforms*;
			python ${params.workpath}/script_concat_rsem.py -o ${params.rsem}/MATRICE_isoforms_FPKM -s FPKM -i ${params.rsem}/*/*isoforms*;
			#Script de concatenation des fichiers de quantifications des genes en matrice globale
			python ${params.workpath}/script_concat_rsem.py -o ${params.rsem}/MATRICE_genes_expeced_count -s expected_count -i ${params.rsem}/*/*genes*;
			python ${params.workpath}/script_concat_rsem.py -o ${params.rsem}/MATRICE_genes_TPM -s TPM -i ${params.rsem}/*/*genes*;
			python ${params.workpath}/script_concat_rsem.py -o ${params.rsem}/MATRICE_genes_FPKM -s FPKM -i ${params.rsem}/*/*genes*;
			echo "\n OK - RSEM sort on Aligned.toTranscriptome.out.bam.ConvertedBam - Run at  `date`" >> ${params.output}/README;
	fi		
    """   
    else 
    """
		echo "RSEM not (finished) yet. \n";
    """  
    }


params.BAMtrans ="${params.workpath}/${params.resultsdir}/STAR_mapping/*/Aligned.toTranscriptome.sort.bam"	
BAMtrans_ch = Channel.fromPath(params.BAMtrans)

Channel
    .fromFilePairs( params.reads )                                             
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }  
    .set { read_pairs3 }

conditions = Channel.value( "${params.conditions}" )
process Cufflinks {

    publishDir "${params.cufflinksDir}"
    
    cpus 8
	executor 'SLURM'

	input:
	file 'params.BAMtrans' from BAMtrans_ch
	set pair_id, file(reads) from read_pairs3

	
    script:
    if( params.mode == 'cufflinks' )
    """
		if [[ "${params.debug}" =~ .*Cuff.*  ]]; then
			echo "\n OK - NO Cufflinks" >> ${params.output}/README;      
		else 	
			cd ${params.workpath}/${params.resultsdir};
			if [ ! -d "${params.cufflinksDir}" ]; then
				mkdir  ${params.cufflinksDir}/; chmod 777  ${params.cufflinksDir}/; 
			fi
			if [ ! -d "${params.cufflinksDir}/${pair_id}/" ]; then
				cd  ${params.cufflinksDir}/; mkdir ${pair_id}/; chmod 777  ${pair_id}/; 
			fi
			cd ${params.cufflinksDir}/${pair_id}/; 
			
			cufflinks --max-intron-length 300000 --num-threads 4 -F 0.1 -j 0.15 --library-type fr-firststrand --GTF-guide ${params.annotation} -o ${params.cufflinksDir}/${pair_id} ${params.mappingDir}/${pair_id}/Aligned.sortedByCoord.out.bam;
			echo "\n OK - Cufflinks on ${params.mappingDir}/${pair_id}/Aligned.sortedByCoord.out.bam - Run at  `date`" >> ${params.output}/README;
			#remove tmp directory
			rm -rf ${params.cufflinksDir}/tmp;
			#quality
			cp ${params.cufflinksDir}/${pair_id}/isoforms.fpkm_tracking ${params.multiQC}/${pair_id}_isoforms.fpkm_tracking;
			cp ${params.cufflinksDir}/${pair_id}/genes.fpkm_tracking ${params.multiQC}/${pair_id}_genes.fpkm_tracking;     
		fi
    """   
    else 
    """
		echo "Cufflinks not (finished) yet.";
    """  
    }


params.GTF ="${params.workpath}/${params.resultsdir}/CUFFLINKS/*/transcripts.gtf"
transcripts_ch = Channel.fromPath(params.GTF)

process Cuffmerge {
    
    cpus 8
	executor 'SLURM'
	

	input:
	file 'params.GTF' from transcripts_ch
	
    script:
    if( params.mode == 'cuffmerge' )
    """
		cd ${params.workpath}/${params.resultsdir};
		#GTF list
		ls -1 ${params.cufflinksDir}/*/transcripts.gtf >> ${params.cufflinksDir}/CufflinksGTF.txt; 
		cd ${params.workpath}/${params.resultsdir};
		if [ ! -d "${params.cuffmergeDir}" ]; then
			mkdir  ${params.cuffmergeDir}/; chmod 777  ${params.cuffmergeDir}/; 
		fi
		cd ${params.cuffmergeDir}/; 
     
		if [[ "${params.debug}" =~ .*Mrg.*  ]]; then
			echo "\n OK - NO Cuffmerge" >> ${params.output}/README;      
		else 	
			#cuffmerge -o ${params.cuffmergeDir}/  -p 4 --ref-sequence  ${params.genome}  -g ${params.annotation}  --min-isoform-fraction 0.05 ${params.cufflinksDir}/CufflinksGTF.txt;
			#cuffmerge -o ${params.cuffmergeDir}/  -p 4 --ref-sequence  ${params.genome}  -g ${params.annotationGFF3}  --min-isoform-fraction 0.05 ${params.cufflinksDir}/CufflinksGTF.txt;
			#without GTF reference (optional)
			cuffmerge -o ${params.cuffmergeDir}/ -p 4 --ref-sequence  ${params.genome} --min-isoform-fraction 0.05 ${params.cufflinksDir}/CufflinksGTF.txt;
			#remove tmp directory
			#rm -rf ${params.cuffmergeDir}/tmp;
			echo "\n OK - Cuffmerge on ${params.cufflinksDir}/CufflinksGTF.txt - Run at  `date`" >> ${params.output}/README;
			#quality
			cp ${params.cuffmergeDir}/logs/run.log ${params.multiQC}/run.log;
		fi	 
    """
    else 
    """
		echo "Cuffmerge not (finished) yet. \n";
    """
    }


params.BAM ="${params.workpath}/${params.resultsdir}/STAR_mapping/*/Aligned.sortedByCoord.out.bam"	
Channel
    .fromFilePairs( params.reads )                                             
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }  
    .set { read_pairs_3 } 

process fCounts {
    
    cpus 8
	executor 'SLURM'

	input:
	file 'params.BAM'
	set pair_id, file(reads) from read_pairs_3

    script:
    if( params.mode == 'fCounts' )
    """
		cd ${params.workpath}/${params.resultsdir};
		if [ ! -d "${params.fCounts}" ]; then
			mkdir  ${params.fCounts}/; chmod 777  ${params.fCounts}/; 
		fi
		cd ${params.fCounts}/; 
		if [[ "${params.debug}" =~ .*FC.*  ]]; then
			echo "\n OK - NO FeatureCounts" >> ${params.output}/README;      
		else 
			featureCounts -s 2 -M -O --primary -t exon -g transcript_id -a ${params.cuffmergeDir}/merged.gtf -o ${params.fCounts}/transcripts.fcounts ${params.mappingDir}/${pair_id}/Aligned.sortedByCoord.out.bam; 
			echo "\n OK - featureCounts on ${params.mappingDir}/${pair_id}/Aligned.sortedByCoord.out.bam - Run at  `date`" >> ${params.output}/README;  
		fi	
    """   
    else 
    """
		echo "featureCounts not (finished) yet. \n";
    """  
    }
    
    
process fCounts_quality {
    
	executor 'SLURM'

    script:
    if( params.mode == 'fCountsQuality' )
    """
		#quality
		cp ${params.fCounts}/transcripts.fcounts.summary ${params.multiQC}/.; 
    """   
    else 
    """
		echo "featureCounts Quality not (finished) yet. \n";
    """  
    }    
    

//params.BAM ="${params.workpath}/${params.resultsdir}/STAR_mapping/*/Aligned.sortedByCoord.out.bam"	
Channel
    .fromFilePairs( params.reads )                                             
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }  
    .set { read_pairs_4 } 

process fCountsOnRef {
    
    cpus 8
	executor 'SLURM'

	input:
	file 'params.BAM'
	set pair_id, file(reads) from read_pairs_4

    script:
    if( params.mode == 'fCountsOnRef' )
    """
		cd ${params.workpath}/${params.resultsdir};
		if [ ! -d "${params.fCountsOnRef}" ]; then
			mkdir  ${params.fCountsOnRef}/; chmod 777  ${params.fCountsOnRef}/; 
		fi
		cd ${params.fCountsOnRef}/; 
		if [[ "${params.debug}" =~ .*fref.*  ]]; then
			echo "\n OK - NO FeatureCounts" >> ${params.output}/README;      
		else 
			featureCounts -s 2 -M -O --primary -t exon -g transcript_id -a ${params.annotation} -o ${params.fCountsOnRef}/transcripts.fcounts ${params.mappingDir}/${pair_id}/Aligned.sortedByCoord.out.bam; 
			echo "\n OK - featureCounts for GTF reference on ${params.mappingDir}/${pair_id}/Aligned.sortedByCoord.out.bam - Run at  `date`" >> ${params.output}/README;
		fi	
    """   
    else 
    """
		echo "featureCounts not (finished) yet. \n";
    """  
    }

process fCountsOnRef_quality {
    
	executor 'SLURM'

    script:
    if( params.mode == 'fCountsQualityOnRef' )
    """
		#quality
		cp ${params.fCountsOnRef}/transcripts.fcounts.summary ${params.multiQC}/.;  
    """   
    else 
    """
		echo "featureCounts Quality not (finished) yet. \n";
    """  
    }    


process FEELnc {
    
    cpus 16
	executor 'SLURM'
	memory '60 GB'
	
	input:
	file 'params.GTF'
	file 'params.fasta'

    script:
    //
    //Step 1 : Feelnc Filter
    //
    if( params.mode == 'FEELnc' )
    """
		cd ${params.workpath}/${params.resultsdir};
		if [ ! -d "${params.feelnc}" ]; then
			mkdir  ${params.feelnc}; chmod 777  ${params.feelnc};
		fi
     
		#step 1 : FILTER
		if [ ! -d "${params.feelnc}/filter" ]; then
			cd  ${params.feelnc}; mkdir filter/; chmod 777 filter/;
		fi
     
		if [[ "${params.debug}" =~ .*FEELnc.*  ]]; then
			echo "\n OK - NO FEELnc filter" >> ${params.output}/README;      
		else 
			FEELnc_filter.pl -i ${params.cuffmergeDir}/merged.gtf -a ${params.annotation}  --biotype transcript_biotype=protein_coding --monoex=1 --size=200  -o ${params.feelnc}/filter/FEELnc_filter.log --proc=12  > ${params.feelnc}/filter/candidate_lncRNA.gtf;
			#quality on filter step
			cp ${params.feelnc}/filter/FEELnc_filter.log 	${params.multiQC}/.;  
		fi
		cd ${params.workpath}/${params.resultsdir};
		#step 2 : PROTEIN PROT shuffle et intergenic
		if [ ! -d "${params.feelnc}/protein_coding/" ]; then
			cd ${params.feelnc}; mkdir protein_coding/; chmod 777 protein_coding/;
		fi
		if [ ! -d "${params.feelnc}/protein_coding/shuffle" ]; then
			cd ${params.feelnc}/protein_coding/ ; mkdir shuffle/; chmod 777 shuffle/;
		fi
		if [ ! -d "${params.feelnc}/protein_coding/intergenic/" ]; then
			cd ${params.feelnc}/protein_coding/; mkdir intergenic/; chmod 777 intergenic/;
		fi
     
		if [[ "${params.debug}" =~ .*FEELnc.*  ]]; then
			echo "\n OK - NO FEELnc codprot" >> ${params.output}/README;      
		else 
			grep "protein_coding" ${params.annotation} > ${params.feelnc}/training_prot.gtf;
			cd ${params.feelnc}/protein_coding/shuffle/; 
			FEELnc_codpot.pl -i ${params.feelnc}/filter/candidate_lncRNA.gtf -a ${params.feelnc}/training_prot.gtf -g ${params.genome} --mode=shuffle;
		 
			#quality on protein coding shuffle step
			cp ${params.feelnc}/protein_coding/shuffle/*/candidate_lncRNA.gtf_RF_learningData.txt 	${params.multiQC}/shuffle_candidate_lncRNA.gtf_RF_learningData.txt; 
			cp ${params.feelnc}/protein_coding/shuffle/*/candidate_lncRNA.gtf_RF_statsLearn_CrossValidation.txt 	${params.multiQC}/shuffle_candidate_lncRNA.gtf_RF_statsLearn_CrossValidation.txt; 
			cp ${params.feelnc}/protein_coding/shuffle/*/candidate_lncRNA.gtf_RF.txt 	${params.multiQC}/shuffle_candidate_lncRNA.gtf_RF.txt; 
			cp ${params.feelnc}/protein_coding/shuffle/*/candidate_lncRNA.gtf_RF_summary.txt 	${params.multiQC}/shuffle_candidate_lncRNA.gtf_RF_summary.txt;      
		fi
	
    cd ${params.feelnc}/protein_coding/intergenic/; 
     
    if [[ "${params.debug}" =~ .*FEELnc.*  ]]; then
		echo "\n OK - NO FEELnc protein coding" >> ${params.output}/README;      
    else 
		FEELnc_codpot.pl -i ${params.feelnc}/filter/candidate_lncRNA.gtf -a ${params.feelnc}/training_prot.gtf -g ${params.genome}  -b transcript_biotype=protein_coding -b transcript_status=KNOWN --mode=intergenic;
		 
		#quality on protein coding intergenic step
		cp ${params.feelnc}/protein_coding/intergenic/*/candidate_lncRNA.gtf_RF_learningData.txt 	${params.multiQC}/intergenic_candidate_lncRNA.gtf_RF_learningData.txt; 
		cp ${params.feelnc}/protein_coding/intergenic/*/candidate_lncRNA.gtf_RF_statsLearn_CrossValidation.txt 	${params.multiQC}/intergenic_candidate_lncRNA.gtf_RF_statsLearn_CrossValidation.txt; 
		cp ${params.feelnc}/protein_coding/intergenic/*/candidate_lncRNA.gtf_RF.txt 	${params.multiQC}/intergenic_candidate_lncRNA.gtf_RF.txt; 
		cp ${params.feelnc}/protein_coding/intergenic/*/candidate_lncRNA.gtf_RF_summary.txt 	${params.multiQC}/intergenic_candidate_lncRNA.gtf_RF_summary.txt;      
    fi


    #step 3 : CLASSIFIER shuffle et intergenic
    if [ ! -d "${params.feelnc}/classifier/" ]; then
		cd ${params.feelnc}/; mkdir classifier/; chmod 777 classifier/;
    fi
    if [ ! -d "${params.feelnc}/classifier/shuffle" ]; then
		cd ${params.feelnc}/classifier/; mkdir shuffle/; chmod 777 shuffle/;
    fi
    if [ ! -d "${params.feelnc}/classifier/intergenic/" ]; then
			cd ${params.feelnc}/classifier/; mkdir intergenic/; chmod 777 intergenic/;
    fi 
     
    if [[ "${params.debug}" =~ .*FEELnc.*  ]]; then
		echo "\n OK - NO FEELnc classifier" >> ${params.output}/README;      
    else 
        #shuffle
        #[Possibly] {INPUT}.noORF.gtf 
		if [ ! -d "${params.feelnc}/protein_coding/shuffle/*/*noORF.gtf" ]; then
			cp ${params.feelnc}/protein_coding/shuffle/*/*lncRNA.gtf  ${params.feelnc}/lst_shuffle_lncRNA_noORF.gtf;
		else
		cat ${params.feelnc}/protein_coding/shuffle/*/*noORF.gtf ${params.feelnc}/protein_coding/shuffle/*/*lncRNA.gtf > ${params.feelnc}/lst_shuffle_lncRNA_noORF.gtf;
		fi
		FEELnc_classifier.pl -i ${params.feelnc}/lst_shuffle_lncRNA_noORF.gtf -a  ${params.feelnc}/training_prot.gtf  > ${params.feelnc}/classifier/shuffle/lncRNA_shuffle_classes.txt;
			
		#intergenic
		#[Possibly] {INPUT}.noORF.gtf  
		if [ ! -d "${params.feelnc}/protein_coding/intergenic/*/*noORF.gtf" ]; then
			cp ${params.feelnc}/protein_coding/intergenic/*/*lncRNA.gtf ${params.feelnc}/lst_intergenic_lncRNA_noORF.gtf;
		else
			cat ${params.feelnc}/protein_coding/intergenic/*/*noORF.gtf ${params.feelnc}/protein_coding/intergenic/*/*lncRNA.gtf > ${params.feelnc}/lst_intergenic_lncRNA_noORF.gtf;
		fi
			 
		FEELnc_classifier.pl -i ${params.feelnc}/lst_intergenic_lncRNA_noORF.gtf -a  ${params.feelnc}/training_prot.gtf   > ${params.feelnc}/classifier/intergenic/lncRNA_intergenic_classes.txt;
			 
		##quality on classifier shuffle and intergenic step
		cp ${params.feelnc}/classifier/*.log ${params.multiQC}/.;  
	fi	
    echo "\n OK - FEELnc (filter, protein coding and classifier) on ${params.cuffmergeDir}/merged.gtf - Run at  `date`" >> ${params.output}/README;
    """
    else 
    """
		echo "Feelnc not (finished) yet. \n";
    """
    }

process quality {
    script:
    //
    //multiQC
    //
    if( params.mode == 'qual' )
    """
	#last Star mapping log
	mv ${params.workpath}/${params.resultsdir}/Log.out ${params.multiQC}/Log.out; 
	cd ${params.multiQC}/; 
	if [[ "${params.debug}" =~ .*Quality.*  ]]; then
		echo "\n OK - NO quality and statistics" >> ${params.output}/README;      
	else 
		sh ${params.workpath}/report_multiqc.sh;
		echo "OK - multiQC report - Run at  `date`" >> ${params.output}/README ;
	fi 
	#Delete tmp directory generated by Nextflow
	rm -rf ${params.workpath}/work/;
    """
    else 
    """
		echo "Sequences list. \n";
    """
    }