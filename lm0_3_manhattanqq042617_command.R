#chroms <- c(23:23)

#for(chrom in chroms){
#   callstring <- paste("nice mosrun -e -q R --slave --args ", chrom, " < doFileConvert2bdoseMafsl072114.R &")
#   print(callstring)
#   system(callstring)
#}

callstring <- paste("Rscript lm0_3_manhattanqq042617.R")
print(callstring)
system(callstring)

