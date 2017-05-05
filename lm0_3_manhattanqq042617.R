# ssh sulgkerberos@chef.compbio.ucsf.edu
# ssh iqint

#################### Top hit Manhattan Plots (NEW 12/08/15 plot pval < 0.00 for Mx) ############################
basedir <- "/netapp/home/sulgkerberos/sulggi/"
setwd(basedir)

# For the output:
# field2 = chromosome
# field3 = rsid
# field5 = position bp on chromosome
# field30 = anova pval

############################To Prepare File to Read In #####################################

#cd data/vault5/sulggi/joel/results/uarto
#awk -F "," '{print $2, $4, $5, $27}' ExomeGWAS_chr1all.bdosemaf.lm0_3mo.lmm > ExomeGWAS_chr1all.bdosemaf.lm0_3mo.lmmsmall
#Then concatenate all the lmm files together 
#cat *.bdosemaf.lm0_3mo.lmmsmall > lm0_3mo.lmm


########################################################################################
# example of one way to make a basic manhattan plot with all chromosomes on one figure.

#tempdat <- read.table("lme4_model3.lmm", sep=" ", as.is=T, header=F)
tempdat <- read.table("joel/results/uarto/lm0_3mo.lmm", sep=" ", as.is=T, header=F)
names(tempdat) <- c("CHR", "SNP", "BP", "P")

head(tempdat)
  
# make  numeric
#P=as.numeric(tempdat$P)
#CHR=as.numeric(tempdat$CHR)
#SNP=tempdat$SNP
#BP=as.numeric(tempdat$BP) 

# sort by CHR and BP
dat <- tempdat[order(tempdat$CHR, tempdat$BP),]
head(dat)
# CHR         SNP    BP          P
# 3455302   1  rs58108140 10583 1.00000000
# 3455303   1 rs180734498 13302 0.13278706
# 3455304   1 rs140337953 30923 0.57872148
# 3455305   1 rs116400033 51479 0.13730468
# 3455306   1 rs150021059 52238 0.06727875
# 3455307   1 rs140052487 54353 0.20943035

# (only plot pval < 0.01)
# manhattangeno <- dat[dat$P < 0.01,]
#keygeno <- dat[dat$P < 0.01,]
keygeno <- tempdat[tempdat$P < 0.00,]
#rm(dat)

dim(keygeno)

head(keygeno)

chroms <- unique(keygeno$CHR)
#pdf("keyHitschrall.pdf")
#for(chrom in chroms){
 # chromDat <- keygeno[keygeno$CHR==chrom,]
 # plot((-1*log10(keygeno[,4])) ~ keygeno[,3], main="All Chromosomes", pch=20, ymin = 2, col=rgb(0,0,1,.3), xlab="position", ylab="-log10(pval) for geno" )
#}
#dev.off()

source("joel/code/qqman.r")

#pdf("Manhattanchrall.pdf")
#manhattan(keygeno)
#dev.off()


######################### QQ PLOTS (Use Genotyped SNPs Only) ##################
## Genotyped SNPs (BIM file) - chromosomes ##
chrallbim <-read.table("joel/code/ExomeGWASmergefinal.bim", header=FALSE)
head(chrallbim)
# V1          V2     V3     V4 V5 V6
# 1  1 rs150690004 0.0007  59359  0  G
# 2  1 rs141776804 0.0007  59453  0  T
# 3  1   rs3094315 0.0000 732292  T  C
# 4  1   rs3131972 0.0000 732446  C  T
# 5  1  rs12562034 0.0000 748174  A  G
# 6  1 rs148989274 0.0084 752348  A  C

genotyped <-  data.frame(chrallbim$V2)
names(genotyped) <-c("SNP")
head(genotyped)
#           SNP
# 1  rs150690004
# 2  rs141776804
# 3    rs3094315

genotyped$TYPED <- rep(1,nrow(genotyped))
head(genotyped)
#         SNP    TYPED
# 1 rs150690004     1
# 2 rs141776804     1
# 3   rs3094315     1
# 4   rs3131972     1
# 5  rs12562034     1
# 6 rs148989274     1

#dim(dat)
dat <- tempdat
dim(dat)
dim(genotyped)

#genodat <- merge(dat,genotyped,by.x="SNP", by.y="SNP", all=TRUE) ## this takes FOREVER!! had to stop it.



library(sqldf)

typedresults <- sqldf(c(
    "create index indSnpDat on dat(SNP)",
    "create index indSnpGenotyped on genotyped(SNP)",
    "select d.* from dat as d join genotyped as g on d.SNP=g.SNP"))

write.csv(typedresults, file="lm0_3mo_typed.lmm", row.names=F)

#typedresults <- genodat[ which(genodat$TYPED==1),]
#write.table(genodat, file=paste(f, ".typed", sep="") , row.names=F, col.names=F, sep=" ", quote=F)

allResults <- sqldf(c(
    "create index indSnpDat on dat(SNP)",
    "create index indSnpGenotyped on genotyped(SNP)",
    "select d.*, TYPED from dat as d left join genotyped as g on d.SNP=g.SNP"))
write.csv(allResults, file="lm0_3mo.lmm", row.names=F)

genodat <- allResults
dim(genodat)

source("joel/code/qqman.r")


#qq(typedresults$P)
#pdf("qqchralltyped.pdf")
#qq(typedresults$P)
#dev.off()

#qq(genodat$P)

#pdf("qqchrall.pdf")
#qq(genodat$P)
#dev.off()

# GENOTYPED SNPS ONLY
typedMedianP <- median(typedresults$P)
typedMedianChi2 <- qchisq(p=typedMedianP, df=1, lower.tail=F)
expectedMedianChi2 <- qchisq(p=.5, df=1, lower.tail=F)
typedLambdaGC <- typedMedianChi2/expectedMedianChi2

# ALL GENOTYPED+IMPUTED SNPS
allMedianP <-	 median(genodat$P)
allMedianChi2	 <- qchisq(p=allMedianP, df=1, lower.tail=F)
#expectedMedianChi2 <- qchisq(p=.5, df=1, lower.tail=F)
allLambdaGC <- allMedianChi2/expectedMedianChi2

# GENOTYPED SNPS QQ PLOT
bitmap("qqtyped_lm0_3mo.png")
qq(typedresults$P)
#pdf("qqtyped_lm0_3mo.pdf")
text(x=1.2, y=1.8, label=sprintf("lambda_gc = %f",round(typedLambdaGC,4) ) )
dev.off()

# ALL GENOTYPED+IMPUTED SNPS QQ PLOT
bitmap("qqtypedimputed_lm0_3mo.png")
qq(genodat$P)
#pdf("qqtypedimputed_lm0_3mo.pdf")
text(x=1.2, y=1.8, label=sprintf("lambda_gc = %f",round(allLambdaGC,4) ) )
dev.off()

# NOW ADJUST GENOTYPED WITH LAMBDAGC
typedPvals <- typedresults$P
typedChisqStatistics <- qchisq(p=typedPvals, df=1, lower.tail=F)
typedAdjustedChisqStatistics <- typedChisqStatistics/typedLambdaGC
typedAdjustedPvals <- pchisq(q=typedAdjustedChisqStatistics, df=1, lower.tail=F)

# ADJUST GENOTYPED WITH LAMBDAGC QQ PLOT
bitmap("qqtypedAdj_chr22_sqrtcd4_regular.png")
qq(typedAdjusted$P)
#pdf("qqtypedAdj_lm0_3mo.pdf")
text(x=1.2, y=1.8, label=sprintf("lambda_gc = %f",round(allLambdaGC,4) ) )
dev.off()


