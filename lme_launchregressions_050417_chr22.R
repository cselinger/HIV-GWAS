#bdosefiles <-dir("/sulggi/newbdoses")[grepl("bdose.$", dir("/sulggi/newbdoses")    )]
bdosefiles <- "ExomeGWAS_chr22all.bdosemaftest"

for(bdosefile in bdosefiles){
    callstring <- paste("Rscript lme_doregressions_050417.R", bdosefile)
    print(callstring)
    system(callstring)
}
