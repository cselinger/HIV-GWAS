#bdosefiles <-dir("/sulggi/newbdoses")[grepl("bdosemaf$", dir("/sulggi/newbdoses")    )]
bdosefiles <- "ExomeGWAS_chr22all.bdosesmaf"

for(bdosefile in bdosefiles){
#   callstring <- paste("R --slave --args ", bdosefile, " < qb3_lmertest_sqrtcd4_doregressions_092816.R &")
    callstring <- paste("Rscript qb3_lmertest_sqrtcd4_doregressions_010517_KRp.R", bdosefile)
    print(callstring)
    system(callstring)
}
