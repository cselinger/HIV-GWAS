bdosefiles <-dir("/netapp/home/sulgkerberos/sulggi/newbdoses")[grepl("bdosemaf$", dir("/netapp/home/sulgkerberos/sulggi/newbdoses")    )]

for(bdosefile in bdosefiles){
  #   callstring <- paste("R --slave --args ", bdosefile, " < qb3_lmertest_sqrtcd4_doregressions_092816.R &")
  callstring <- paste("Rscript lm0_3mo_doregressions042617.R", bdosefile)
  print(callstring)
  system(callstring)
}