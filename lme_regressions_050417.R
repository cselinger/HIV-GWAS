 
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

# lmer gives "t-values", but df not obviously defined, so for now pretend they are z-values to get p-values
z2p <- function(z){
  p <- 2*pnorm((-1)*abs(z) )
  return(p)
}

z2pv <- function(zv){
  ret <- sapply(zv, z2p)
  return(ret)
}

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

     #   i = 1
     for(i in 1:nr){
       fields <- infields[[i]]              
       bdosemafHeader <- fields[1:lengthbdoseheader]
       # lenth of results: 3 + 1 + (number genetic terms)*5
       ngeneticTerms <- 4  # modify for other models- geno, timeint90:geno, timeint90365:geno, timeint365plus:geno
       results <- rep("", 3+1+5*ngeneticTerms)
       
       try({
         bdoses <- as.numeric(fields[(lengthbdoseheader+1):length(fields)] )
         geno <- bdoses[phenoUse$genoOrder]
         regData <- cbind(phenoUse, geno=geno, stringsAsFactors=F)
         
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

            lmod <- as.numeric(logLik(mod) )
            #coefs <- summary(mod, type=2, ddf="Kenward-Roger")$coef
            coefs <- summary(mod)$coef
            
            coefsSmall <- coefs[grepl("geno", row.names(coefs) ) , ]
            coefsSmall <- as.matrix(coefsSmall)    
            
         # new
         # if there is only 1 genetic term, the previous line returns a vector of coefficients without row names
         # rather than a matrix with row names
         if(( class(coefsSmall) != "matrix") && (ngeneticTerms == 1) ){
           coefsSmall <- t(as.matrix(coefsSmall))
           row.names(coefsSmall) <- c("geno")
         }
            print(coefsSmall)
            # end new

         pvalsSmall <- z2pv(coefsSmall[,3])  # convert t-values to p-values (normal approx)
       
         basemod <- lmerTest::lmer(sqrtcd4art ~  (1 | id) + timeint90 + timeint90365 + timeint365plus 
                     + age + timeint90:age + timeint90365:age + timeint365plus:age 
                     + factor(sex) + timeint90:factor(sex) + timeint90365:factor(sex) + timeint365plus:factor(sex)
                     + sqrtcd4base + timeint90:sqrtcd4base + timeint90365:sqrtcd4base + timeint365plus:sqrtcd4base
                     + pca1 + pca2 + pca3 + pca4 + pca5 
                     + pca6 + pca7 + pca8 + pca9 + pca10,
                     na.action=na.omit, REML=F,
                     data=regData[!incomplete & regData$cohort=="UARTO",] );
         
        lbasemod <- as.numeric(logLik(basemod) )
#        anovapval <- anova(basemod, mod)[2,7]

results <- c(lmod, lbasemod, nrow(coefsSmall), row.names(coefsSmall), t(coefsSmall))
       })
       resultsLine <- c(bdosemafHeader, results)
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
}
