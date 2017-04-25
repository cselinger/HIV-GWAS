##make sure that you have all packages needed installed
packages <- c("ggplot2", "data.table", "lmerTest","splines","pbkrtest","optimx")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}

#load packages
invisible(lapply(packages, library, character.only=T,quietly=T))



setwd('c:/users/cselinger/Dropbox (IDM)/HIV-GWAS/')#add your dropbox path....

load('Sulggi/code/cd4Setupsleigen101116.rdata')
phenoUse<-as.data.table(phenoUse)
phenoUse<-phenoUse[cohort=="UARTO",]
phenoUse<-phenoUse[,sqrtCD4:=sqrt(cd4)]
phenoUse<-phenoUse[!id%in% unique(phenoUse[sqrtCD4>35,id]),]#keep only ids with sqrt cd4 <35
phenoUse$sex<-as.factor(phenoUse$gender);phenoUse[,gender:=NULL];levels(phenoUse$sex)<-c('Male','Female')
phenoUse$id<-as.factor(phenoUse$id)
phenoUse$age<-cut(phenoUse$age,c(seq(15,50,10),100))




blocksize=1
lengthbdoseheader <- 21
bdosefile <- "ExomeGWAS_chr22all.bdosesmaf"
infilename <- paste0('Sulggi/newbdoses/',bdosefile)
infile <- file(infilename, open="r")
inblock <- readLines(con=infile, n=blocksize)
infields <- sapply(inblock, strsplit, " ")
nr <- length(inblock)


i=nr
fields <- infields[[i]]              
bdosemafHeader <- fields[1:lengthbdoseheader]
ngeneticTerms <- 4  # modify for other models- geno, treatdays:geno, treatdayspp90:geno, treatdayspp365:geno
results <- rep("", 3+1+5*ngeneticTerms)

bdoses <- as.numeric(fields[(lengthbdoseheader+1):length(fields)] )
geno <- bdoses[phenoUse$genoOrder]
regData <- cbind(phenoUse, geno=geno, stringsAsFactors=F)
regData$geno <- as.character(regData$geno)


regData<-as.data.frame(regData)
regData$knot<- knots(regData$treatdays, c(0,90,365))
regData<-as.data.table(regData)

incomplete <- with(regData,
                   (is.na(sqrtCD4) | is.na(treatdays) | is.na(id) | is.na(geno) |
                      is.na(sex) | is.na(cohort) |
                      is.na(treatdayspp90) | is.na(treatdayspp365) |
                      is.na(timeint90) | is.na(timeint90365) | is.na(timeint365plus) |
                      is.na(pca1)     | is.na(pca2) | is.na(pca3) |
                      is.na(pca4)     | is.na(pca5) | is.na(pca6) |
                      is.na(pca7)     | is.na(pca8) | is.na(pca9) |
                      is.na(pca10)    | is.na(genoOrder) |
                      is.na(knot.0)   | is.na(knot.90)  | is.na(knot.365)) )




m1 <- lm(sqrtCD4~ age+sex+knot, data=regData[!incomplete,] )
matX<-expand.grid(treatdays=c(0,90,365,ceiling(max(regData$treatdays))),age=factor(levels(regData$age)[1:4]),sex=factor(levels(regData$sex)))
matX$knot <- knots(matX$treatdays, c(0,90,365))
matX$predict <- predict(m1, newdata=matX)



# plot of data and predicted values
plot1 <- xyplot(sqrtCD4~treatdays, group=id, data=regData[!incomplete,], type="b")
plot2 <- xyplot(predict~treatdays, data=matX, type="b", col="black", lwd=3)
plot1+plot2


ggplot(matX)+geom_line(aes(x=treatdays,y=predict))+facet_grid(age~sex)



m1 <- lm(sqrtCD4~ age+geno+knot, data=regData[!incomplete,] )
matX<-expand.grid(treatdays=c(0,90,365,ceiling(max(regData$treatdays))),age=factor(levels(regData$age)[1:4]),geno=factor(unique(regData$geno)))
matX$knot <- knots(matX$treatdays, c(0,90,365))
matX$predict1 <- predict(m1, newdata=matX)

m2 <- lmer(sqrtCD4~ (1 | id) + age+geno+knot, data=regData[!incomplete,])
matX$predict2 <- predict(m2, re.form=NA,newdata=matX)

ggplot(matX)+geom_line(aes(x=treatdays,y=predict^2))+ggtitle('Fixed effect model predictions')+xlab('days on treatment')+ylab("CD4 T cell counts")+facet_grid(age~geno)
ggplot(regData[!incomplete,] )+geom_line(aes(x=treatdays,y=sqrtCD4,color=id))+facet_grid(age~geno)+scale_color_discrete(guide = FALSE)




matXX<-merge(regData[!incomplete,],matX,by=c('age','geno'))
ggplot(matXX)+geom_line(aes(x=treatdays.x,y=sqrtCD4,color=id))+facet_grid(age~geno)+scale_color_discrete(guide = FALSE)+
  geom_line(aes(x=treatdays.y,y=predict1))

ggplot(matXX)+geom_line(aes(x=treatdays.x,y=sqrtCD4,color=id))+facet_grid(age~geno)+scale_color_discrete(guide = FALSE)+
  geom_line(aes(x=treatdays.y,y=predict2))




m2 <- lmer(sqrtCD4~ (1 | id) + age+geno+knot, data=regData[!incomplete,])
matX$predict2 <- predict(m2, re.form=NA,newdata=matX)




m2 <- lmer(sqrtCD4~ (1 + knot | id) + age+geno, data=regData[!incomplete,])
m2 <- lmer(sqrtCD4~ (1|id) + (knot -1 | id) + age+geno, data=regData[!incomplete,])


coefs <- coef(summary(m2))
coefsSmall <- coefs[grepl("geno", row.names(coefs) ) , ]
##kenward-roger
require(pbkrtest)
df.KR<-get_ddf_Lb(m2,fixef(m2))
coefsSmall<-as.data.frame(coefsSmall)
coefsSmall$df.KR<-df.KR
coefsSmall$pvalsKR<-2*(1-pt(abs(coefsSmall[,'t value']),df.KR))        





mod <- lmer(sqrtCD4 ~  (1 | id) + treatdays
            + treatdayspp90 + treatdayspp365
            + geno +  timeint90:geno 
            + timeint90365:geno + timeint365plus:geno
            + sex + pca1 + pca2 + pca3 + pca4 + pca5 
            + pca6 + pca7 + pca8 + pca9 + pca10,
            na.action=na.omit, REML=F,
            data=regData[!incomplete & regData$cohort=="UARTO",] );



  