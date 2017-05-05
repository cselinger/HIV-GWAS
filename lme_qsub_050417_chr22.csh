#!/bin/bash                         
#$ -S /bin/bash
#$ -o /netapp/home/sulgkerberos/sulggi/joel/results/uarto
#$ -e /netapp/home/sulgkerberos/sulggi/joel/results/uarto
#$ -cwd                            
#$ -r y                            
#$ -j y                            
#$ -l mem_free=10G                  
#$ -l arch=linux-x64              
#$ -l netapp=10G,scratch=10G         
#$ -l h_rt=336:00:00  
#$ -t 1

cd /netapp/home/sulgkerberos/sulggi/joel/code

name -a
date
echo "Job ID: $SGE_JOBID"
echo "Task ID: $SGE_TASK_ID"

# Read task parameters (this is where list all the chr to run)
# and define them as environmental variables PARAM1-2 (PARAM are COLUMNS in lme_param_chr22.txt)
# above #$ -t 1 (refers to the LINES in lme_param_chr22.txt)

export $(head -$SGE_TASK_ID lme_param_chr22.txt | tail -1 | awk 'BEGIN{ OFS=""; } {printf "PARAM1="; print $1; }')
                               
Rscript $PARAM1 $PARAM2 "/netapp/home/sulgkerberos/sulggi/joel/results/uarto/" "/netapp/home/sulgkerberos/sulggi/newbdoses/" 


