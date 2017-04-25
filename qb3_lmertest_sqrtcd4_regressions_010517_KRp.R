 
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


lme4cd4sqrtregressions <- function(phenoUse, indir, genofileRootName, outdir, blocksize=25, counters=T){
  infilename <- paste(indir, genofileRootName, sep="")
  outfilename <- paste(outdir, genofileRootName, ".qb3_lmertest_sqrtcd4_KRp.lmm", sep="")
  
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
       ngeneticTerms <- 4  # modify for other models- geno, treatdays:geno, treatdayspp90:geno, treatdayspp365:geno
       results <- rep("", 3+1+5*ngeneticTerms)
       
       try({
         bdoses <- as.numeric(fields[(lengthbdoseheader+1):length(fields)] )
         geno <- bdoses[phenoUse$genoOrder]
         regData <- cbind(phenoUse, geno=geno, stringsAsFactors=F)
         
       #regData <- regData[regData$monthart <= maxTreatMonths,]
       sex <- rep(NA, nrow(regData))
       sex[regData$gender=="1:Male"] <- 0
       sex[regData$gender=="2:Female"] <- 1
       regData$sex <- sex
       #names(regData)[9] <- "cd4" in my dataset, column is 7, not 9, and it is already called cd4
       regData$sqrtCD4 <- sqrt(regData$cd4)
       
       incomplete <- with(regData,
                          (is.na(sqrtCD4) | is.na(treatdays) | is.na(id) | is.na(geno) |
                               is.na(sex) | is.na(cohort) |
                               is.na(treatdayspp90) | is.na(treatdayspp365) |
                               is.na(timeint90) | is.na(timeint90365) | is.na(timeint365plus) |
                               is.na(pca1)     | is.na(pca2) | is.na(pca3) |
                               is.na(pca4)     | is.na(pca5) | is.na(pca6) |
                               is.na(pca7)     | is.na(pca8) | is.na(pca9) |
                               is.na(pca10)    | is.na(genoOrder) ) )

mod <- lmer(sqrtCD4 ~  (1 | id) + treatdays
            + treatdayspp90 + treatdayspp365
            + geno +  timeint90:geno 
            + timeint90365:geno + timeint365plus:geno
            + sex + pca1 + pca2 + pca3 + pca4 + pca5 
            + pca6 + pca7 + pca8 + pca9 + pca10,
            na.action=na.omit, REML=F,
            data=regData[!incomplete & regData$cohort=="UARTO",] );

            lmod <- as.numeric(logLik(mod) )
            #coefs <- summary(mod, type=2, ddf="Kenward-Roger")$coef
            coefs <- coef(summary(mod))
            
            coefsSmall <- coefs[grepl("geno", row.names(coefs) ) , ]
            coefsSmall <- as.matrix(coefsSmall)    
            
         # new
         # if there is only 1 genetic term, the previous line returns a vector of coefficients without row names
         # rather than a matrix with row names
         if(( class(coefsSmall) != "matrix") && (ngeneticTerms == 1) ){
           coefsSmall <- t(as.matrix(coefsSmall))
           row.names(coefsSmall) <- c("geno")
         }
            #print(coefsSmall)
            # end new
      
         #pvalsSmall <- z2pv(coefsSmall[,4])  # convert t-values to p-values (normal approx)
         #pvalsSmall<- coefsSmall[,5]#these are satterthwaite approximation p values, lmer used here is from lmerTest (not form lme4)

         ##kenward-roger
         require(pbkrtest)
         df.KR<-get_ddf_Lb(mod,fixef(mod))
         coefsSmall<-as.data.frame(coefsSmall)
         coefsSmall$df.KR<-df.KR
         coefsSmall$pvalsKR<-2*(1-pt(abs(coefsSmall[,'t value']),df.KR))        


         #now go back to lmerTest
         mod <- lmerTest::lmer(sqrtCD4 ~  (1 | id) + treatdays
                               + treatdayspp90 + treatdayspp365
                               + geno +  timeint90:geno 
                               + timeint90365:geno + timeint365plus:geno
                               + sex + pca1 + pca2 + pca3 + pca4 + pca5 
                               + pca6 + pca7 + pca8 + pca9 + pca10,
                               na.action=na.omit, REML=F,
                               data=regData[!incomplete & regData$cohort=="UARTO",]
                               #,control = lmerControl(calc.derivs = FALSE) ##this could speed up optimization but does not check for convergence.
         );
         
         coefs <- summary(mod)$coef
         
         coefsSmall1 <- coefs[grepl("geno", row.names(coefs) ) , ]
         coefsSmall1 <- as.matrix(coefsSmall1)
         coefsSmall$pvals.S<-coefsSmall1[,5]
         
         
         print(coefsSmall)
basemod <- lmer(sqrtCD4 ~  geno + (1 | id) + treatdays
                + treatdayspp90 + treatdayspp365 
                + sex + pca1 + pca2 + pca3 + pca4 
                + pca5 + pca6 + pca7 + pca8 + pca9 + pca10,
                na.action=na.omit, REML=F,
                data=regData[!incomplete & regData$cohort=="UARTO",] );

        lbasemod <- as.numeric(logLik(basemod) )
#        anovapval <- anova(basemod, mod)[2,7]
#        anovapval <- anova(basemod, mod, ddf = "Kenward-Roger")[2,8]
       
# results <- c(lmod, lbasemod,anovapval, nrow(coefsSmall), row.names(coefsSmall), t(coefsSmall), pvalsSmall)
#     })
       
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
