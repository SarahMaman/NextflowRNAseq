#!/bin/bash
#SBATCH -p workq
#SBATCH -t 10 #Acceptable time formats include "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds". 

#Load modules
module load bioinfo/multiqc-v1.5


multiqc .


