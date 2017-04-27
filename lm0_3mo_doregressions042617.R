basedir <- "/netapp/home/sulgkerberos/sulggi/"
setwd(basedir)
load( "joel/code/cd4Setupsleigen042617.rdata")
source("joel/code/lm0_3mo_regressions042617.R") # get regressionsLmmUartoKT() function

args <- commandArgs(T)
bdosefile <- args[1]
indir <- paste(basedir, "newbdoses/", sep="")
outdir <- paste(basedir, "joel/results/uarto/", sep="")
#rootname <- paste("ExomeGWAS_chr", chrom, "all", sep="")
 
lmregressions(phenoUse=phenoUse, indir=indir, genofileRootName=bdosefile, outdir=outdir, blocksize=10000, counters=T)
