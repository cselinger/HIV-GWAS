basedir <- "/netapp/home/sulgkerberos/sulggi/"
setwd(basedir)
load( "joel/code/cd4Setupsleigen050417.rdata")
source("joel/code/lme_regressions_050417.R") # get regressionsLmmUartoKT() function

args <- commandArgs(T)
bdosefile <- args[1]
indir <- paste(basedir, "newbdoses/", sep="")
outdir <- paste(basedir, "joel/results/uarto/", sep="")
#rootname <- paste("ExomeGWAS_chr", chrom, "all", sep="")
 
lmeregressions(phenoUse=phenoUse, indir=indir, genofileRootName=bdosefile, outdir=outdir, blocksize=1000, counters=T)
