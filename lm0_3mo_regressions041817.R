
## this is set up to run the univariate lmm for logktrationmum
## in the function regressionsLmmSulggi(), 
## update the following to run other models:
## outfilename
## mod
## basemod
## ngeneticTerms

## phenotype file and correspondence with genotypes prepared as follows: see /Users/sulggi/Sulggi/Imputation/code/phenoUse060714.R
# pheno <- read.csv("/Users/sulggi/Sulggi/Data/UARTOfinalktratio060714eigen.csv", as.is=T)
# samplefile <- read.table("/Users/sulggi/Sulggi/Imputation/chromogenofiles/ExomeGWASmergefinal_chr20.sample", as.is=T, header=T)
# samplefile <- samplefile[2:nrow(samplefile),]
#so bdoseMaf line [1:11]:  mafAboveThreshold, chr, snpid, rsid, position, alleleA, alleleB, maf, baf, minAllele, minAlleleValue, [then bdoses]
# genoOrder <- data.frame(id1=samplefile$ID_1, id2=samplefile$ID_2, genoOrder=1:nrow(samplefile) )
# phenox <- pheno[!is.na(pheno$EXOMEKT),]
# library(sqldf)
# phenoxo <- sqldf("select phenox.*, genoOrder from phenox left join genoOrder on phenox.id=genoOrder.id2")
# phenoUse <- phenoxo[, c("id", "genoOrder", "logktrationmum", "logktrbasenmum", "logktr6nmum", "logktr12nmum", "logktr612nmum", "month", "monthcat", "monthart", "monthcatart", "cohortnum", "gendernum", "EXOMEKT", "pregnant", "pregmen", "logage", "sqrtcd4base", "logvlbase", "pca1", "pca2", "pca3", "pca4", "pca5", "pca6", "pca7", "pca8", "pca9", "pca10")]
# save(pheno, phenox, phenoxo, phenoUse, samplefile, file="ktratioSetupsleigen.rdata")

#library(lme4)
options(scipen=20)

# lmer gives "t-values", but df not obviously defined, so for now pretend they are z-values to get p-values
z2p <- function(z){
  p <- 2*pnorm((-1)*abs(z) )
  return(p)
}

z2pv <- function(zv){
  ret <- sapply(zv, z2p)
  return(ret)
}

## codedir <- "/Users/joel/Desktop/sulggiSept/sept30/newcode/"
## indir <- "/Users/joel/Desktop/sulggiSept/sept30/fromLodz/"
## outdir <- "/Users/joel/Desktop/sulggiSept/sept30/newcode/regressionOutput/"
## genofileRootName <- "test10"
## blocksize <- 1000
## counters <- T
#load( sprintf("%sktratioSetupsleigen.rdata", codedir))

## # each line: (model 1 -> 3 genetic terms)
## 11 fields from bdose file, then:
## number of genetic terms
## names of genetic terms
## 3 fields for each term (estimate, sterr, tvalue)for term1, (estimate, sterr, tvalue) for term2 ...
## pvalues for each of the genetic terms (calculated from t-values, normal approximation)

lmregressions <- function(phenoUse, indir, genofileRootName, outdir, blocksize=1000, counters=T){
  infilename <- paste(indir, genofileRootName, sep="")
  # outfilename <- paste(outdir, genofileRootName, ".model1.lmm", sep="")
  # outfilename <- paste(outdir, genofileRootName, "_model2.lmm", sep="")
   outfilename <- paste(outdir, genofileRootName, ".lm0_3mo.lmm", sep="")
  
  infile <- file(infilename, open="r")
  outfile <- file(outfilename, open="w")
  done <- F
  blocki <- 0
  lengthbdoseheader <- 21
 
  
  while(!done){
     outlines <- c()
     blocki <- blocki + 1
     inblock <- readLines(con=infile, n=blocksize)
     nr <- length(inblock)
     if(nr < blocksize){ done <- T;}
     if(nr == 0){done <- T; break()}
     infields <- sapply(inblock, strsplit, " ")

     #   i = 1
     for(i in 1:nr){
       fields <- infields[[i]]              
       bdosemafHeader <- fields[1:lengthbdoseheader]
       # lenth of results: 3 + 1 + (number genetic terms)*5
       ngeneticTerms <- 1  # modify for other models
       results <- rep("", 3+1+5*ngeneticTerms)
       
       try({
         bdoses <- as.numeric(fields[(lengthbdoseheader+1):length(fields)] )
         geno <- bdoses[phenoUse$genoOrder]
         regData <- cbind(phenoUse, geno=geno, stringsAsFactors=F)
         
#         incomplete <- with(regData, (is.na(logktr6nmum) | is.na(pregmen) | is.na(cohortnum) | is.na(monthart) |
#         is.na(geno) | is.na(pca1) | is.na(pca2) | is.na(pca3) | is.na(pca4) | is.na(pca5) | is.na(pca6) | 
#         is.na(pca7) | is.na(pca8) | is.na(pca9) | is.na(pca10) | is.na(id)))

          incomplete <- with(regData, (is.na(delta_0_3uart) | is.na(sqrtcd4_0uart)| is.na(geno) |
          is.na(age) | is.na(sex) |  is.na(cohort)| is.na(pca1) | is.na(pca2) | is.na(pca3) | is.na(pca4) | 
          is.na(pca5) | is.na(pca6) | is.na(pca7) | is.na(pca8) | is.na(pca9) |is.na(pca10) ) )         
            
                  
#         mod <- lm(logktr6nmum ~  factor(pregmen) + factor(cohortnum) + geno + 
#         pca1 + pca2 + pca3 + pca4 + pca5 + pca6 + pca7 + pca8 + pca9 + pca10, 
#         na.action=na.omit, data=regData[!incomplete,]);
                  
          mod <- lm(delta_0_3uart ~  sqrtcd4_0uart + geno  + sex + age 
                 + pca1 + pca2 + pca3 + pca4 + pca5 
                 + pca6 + pca7 + pca8 + pca9 + pca10,
                 na.action=na.omit, REML=F,
                 data=regData[!incomplete & regData$cohort=="UARTO",] );
                  
         lmod <- as.numeric(logLik(mod) )
         coefs <- summary(mod)$coef
         coefsSmall <- coefs[grep("geno", row.names(coefs) ) , ]
         
         # new
         # if there is only 1 genetic term, the previous line returns a vector of coefficients without row names
         # rather than a matrix with row names
         if(( class(coefsSmall) != "matrix") && (ngeneticTerms == 1) ){
           coefsSmall <- t(as.matrix(coefsSmall))
           row.names(coefsSmall) <- c("geno")
         }
         # end new
# end new
        #pvalsSmall <- z2pv(coefsSmall[,3])  # convert t-values to p-values (normal approx)
        pvalsSmall <- coefsSmall[,4]  # lm() should already have this for your
        coefsSmall <- coefsSmall[,1:3]  # because we were using lmer that gave 
                                        # 3 col coefs matrices, but lm() gives
                                        # and we want to keep code and output
                                        # familiar ...4, and             
        # baseModel: no genotype terms
        #basemod <- lmer(logktrationmum ~  factor(pregmen) + factor(cohortnum) + monthart + 
        #pca1 + pca2 + pca3 + pca4 + pca5 + pca6 + pca7 + pca8 + pca9 + pca10 + (1|id), 
        #na.action=na.omit, data=regData[!incomplete,] , REML=F);

#        basemod <- lm(logktr6nmum ~ factor(pregmen) + factor(cohortnum)
#        + pca1 + pca2 + pca3 + pca4 + pca5 + pca6 + pca7 + pca8 + pca9 +
#        pca10, na.action=na.omit, data=regData[!incomplete,] );

        basemod <- lm(delta_0_3uart ~  sqrtcd4_0uart  + sex + age
                      + pca1 + pca2 + pca3 + pca4 + pca5 
                      + pca6 + pca7 + pca8 + pca9 + pca10,
                      na.action=na.omit, REML=F,
                      data=regData[!incomplete & regData$cohort=="UARTO",] );
        
        lbasemod <- as.numeric(logLik(basemod) )
        anovapval <- anova(basemod, mod)[2,7]
        results <- c(lmod, lbasemod,anovapval, nrow(coefsSmall), row.names(coefsSmall), t(coefsSmall), pvalsSmall)
     })

       resultsLine <- c(bdosemafHeader, results)
       outline <- paste(resultsLine, sep="", collapse=",")
       outlines <- c(outlines,outline )
      
     }
     
     writeLines(outlines, con=outfile)
     if(counters){
       print(c(infilename, blocki*blocksize))
     }
   }
  close(infile)
  close(outfile)
}



