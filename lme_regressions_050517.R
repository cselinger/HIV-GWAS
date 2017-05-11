 
## update the following to run other models:
## outfilename
## mod
## basemod
## ngeneticTerms


library(lme4)
options(scipen=100)
library(lmerTest)
options(ddf="lme4")
library(pbkrtest)
library(data.table)

KenwardRoger=F#if F, don't use KR approximation of df


lmeregressions <- function(phenoUse, indir, genofileRootName, outdir, blocksize=1000, counters=T){
  infilename <- paste(indir, genofileRootName, sep="")
  outfilename <- paste(outdir, genofileRootName, ".lme_regressions.lmm", sep="")
  
  infile <- file(infilename, open="r")
  outfile <- file(outfilename, open="w")

  #maxTreatDays <- 365   # use this if want to only want to look at 1 year and after for regression as per Peter's suggestion
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

     i = 1
     for(i in 1:nr){
       fields <- infields[[i]]              
       bdosemafHeader <- fields[1:lengthbdoseheader]
       # lenth of results: 3 + 1 + (number genetic terms)*5
       ngeneticTerms <- 4  # modify for other models- geno, timeint90:geno, timeint90365:geno, timeint365plus:geno
       results <- rep("", 3+1+5*ngeneticTerms)
       
       try({
         bdoses <- as.numeric(fields[(lengthbdoseheader+1):length(fields)] ) ###should take baf column, should be column #19
         geno <- bdoses[phenoUse$genoOrder]
         regData <- cbind(phenoUse, geno=geno, stringsAsFactors=F)
         
         regData<-data.table(regData)##need to install and load package data.table
         regData[,sqrtcd4base:=na.omit(sqrtcd4base),by='id']#otherwise regData will be empty when removing incomplete rows
         
       #regData <- regData[regData$treatdays >= maxTreatDays,]  # use this if want to only want to look at 1 year and after for regression as per Peter's suggestion 

       incomplete <- with(regData,
                          (is.na(sqrtcd4art) | is.na(id) | is.na(geno) | is.na(timeint90) | is.na(timeint90365) | is.na(timeint365plus) | 
                             is.na(age) | is.na(sex) | is.na(cohort) | is.na(sqrtcd4base) | 
                               is.na(pca1)     | is.na(pca2) | is.na(pca3) |
                               is.na(pca4)     | is.na(pca5) | is.na(pca6) |
                               is.na(pca7)     | is.na(pca8) | is.na(pca9) |
                               is.na(pca10)    | is.na(genoOrder) ) )

       mod <- lmerTest::lmer(sqrtcd4art ~  (1 | id) + geno + timeint90 + timeint90365 + timeint365plus 
                   + timeint90:geno + timeint90365:geno + timeint365plus:geno
                   + age + timeint90:age + timeint90365:age + timeint365plus:age 
                   + factor(sex) + timeint90:factor(sex) + timeint90365:factor(sex) + timeint365plus:factor(sex)
                   + sqrtcd4base + timeint90:sqrtcd4base + timeint90365:sqrtcd4base + timeint365plus:sqrtcd4base
                   + pca1 + pca2 + pca3 + pca4 + pca5 
                   + pca6 + pca7 + pca8 + pca9 + pca10,
                   na.action=na.omit, REML=F,
                   data=regData[!incomplete & regData$cohort=="UARTO",] );

       coefs <- summary(mod)$coef
       coefs <- coefs[grepl("geno", rownames(coefs) ) , ]; foo<-rownames(coefs)
       coefs <- data.table(coefs);coefs$covariate=foo
       colnames(coefs)[grep('Pr(>|t|)',colnames(coefs))] <- 'pvalSW'
       if (KenwardRoger==T){
         coefs$dfKR<-get_ddf_Lb(mod,fixef(mod))
         coefs[,pvalKR:=2*(1-pt(abs(get('t value')),dfKR))]
       }

            
         # new
         # if there is only 1 genetic term, the previous line returns a vector of coefficients without row names
         # rather than a matrix with row names
         if(( class(coefs) != "matrix") && (ngeneticTerms == 1) ){
           coefs <- t(as.matrix(coefs))
           row.names(coefs) <- c("geno")
         }
            print(coefs)
            # end new

       })
       ##write CHR SNP BP PvaluesOfGeno
       resultsLine <- c(bdosemafHeader[c(3,4,5)], as.numeric(unlist(coefs[,grep('pval',colnames(coefs)),with=F])))
       
       ##here are the column names of the output for QQ and manhattan plot
       NAMES<-NULL
       for (j in grep('pval',colnames(coefs))){
         NAMES<-c(NAMES,apply(cbind(coefs$covariate,colnames(coefs)[j]), 1, paste, collapse="_"))
       }
       NAMES<-c('CHR','SNP','BP',NAMES)
       
       
       
       outline <- paste(resultsLine, sep="", collapse=",")
       outlines <- c(outlines,outline )
       
     } # End of SNP 
     
     writeLines(outlines, con=outfile)
     if(counters){
       print(c(infilename, blocki*blocksize))
     }
  } # End of block (not done) 
  
  close(infile)
  close(outfile)
  
  save(NAMES,file=paste(outdir, genofileRootName, ".lme_regressions.NAMES.Rdata", sep=""))
}
