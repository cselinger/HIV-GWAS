# Create PhenoUse for CD4 Recovery Analysis

#install.packages("sqldf")

pheno <- read.csv("/Users/sulggi/Sulggi/GWAS_CD4_Recovery/data/modified/pheno_cd4_042617newpc.csv", as.is=T)
samplefile <- read.table("/Users/sulggi/Sulggi/Imputation/chromogenofiles/ExomeGWASmergefinal_chr20.sample", as.is=T, header=T)
samplefile <- samplefile[2:nrow(samplefile),]
#so bdoseMaf line [1:11]:  mafAboveThreshold, chr, snpid, rsid, position, alleleA, alleleB, maf, baf, minAllele, minAlleleValue, [then bdoses]
genoOrder <- data.frame(id1=samplefile$ID_1, id2=samplefile$ID_2, genoOrder=1:nrow(samplefile) )
phenox <- pheno[!is.na(pheno$EXOMECD4),]
library(sqldf)
phenoxo <- sqldf("select phenox.*, genoOrder from phenox left join genoOrder on phenox.id=genoOrder.id2")
phenoUse <- phenoxo[, c("id", "genoOrder", "cohort", "EXOMECD4", "arv_start_date", "sex", "cd4", "treatdays", 
                        "treatdayspp90", "treatdayspp365", "timeint90", "timeint90365", "timeint365plus", 
                        "pca1", "pca2", "pca3", "pca4", "pca5", "pca6", "pca7", "pca8", "pca9", "pca10",  
                        "age", "preart_cd4", "preart_vl", "vl_result", "cd4_date",
                        "month", "monthcat", "sqrtcd4", "sqrtcd4_0uart", "sqrtcd4_3uart", "sqrtcd4_6uart", 
                        "sqrtcd4_9uart", "sqrtcd4_12uart", "sqrtcd4_24uart","sqrtcd4_maxuart", 
                        "delta_0_3uart", "delta_3_12uart", "delta_12_24uart", "delta_24_maxuart",
                        "sqrtcd4_0arks", "sqrtcd4_3arks", "sqrtcd4_6arks", "sqrtcd4_9arks", "sqrtcd4_12arks", 
                        "delta_0_3arks", "delta_3_12arks", "first", "age_base", "sex_base", 
                        "pca1_base", "pca2_base", "pca3_base", "pca4_base", "pca5_base", 
                        "pca6_base", "pca7_base", "pca8_base", "pca9_base", "pca10_base")]
save(pheno, phenox, phenoxo, phenoUse, samplefile, file="cd4Setupsleigen042617.rdata")