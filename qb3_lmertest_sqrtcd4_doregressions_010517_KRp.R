#basedir <- "/netapp/home/sulgkerberos/sulggi/"
# basedir<- '/home/cselinger/2ndDrive/HIV-GWAS/Sulggi/'

basedir<-ifelse(Sys.info()['user']!='cselinger',"/netapp/home/sulgkerberos/sulggi/",'c:/users/cselinger/Dropbox (IDM)/HIV-GWAS/Sulggi/')

setwd(basedir)

#load( "joel/code/cd4Setupsleigen101116.rdata")
#source("joel/code/qb3_lmertest_sqrtcd4_regressions_010517_KRp.R") # get regressionsLmmUartoKT() function

load( "code/cd4Setupsleigen101116.rdata")
source("code/qb3_lmertest_sqrtcd4_regressions_010517_KRp.R") # get regressionsLmmUartoKT() function


args <- commandArgs(T)
bdosefile <- args[1]
indir <- paste(basedir, "newbdoses/", sep="")
outdir <- paste(basedir,"results/uarto/",sep="")

#outdir <- paste(basedir, "joel/results/uarto/", sep="")
#rootname <- paste("ExomeGWAS_chr", chrom, "all", sep="")
 
lme4cd4sqrtregressions(phenoUse=phenoUse, indir=indir, genofileRootName=bdosefile, outdir=outdir, blocksize=25, counters=T)
