##make sure that you have all packages needed installed
packages <- c("ggplot2", "data.table", "lmerTest","splines","pbkrtest","optimx")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}

#load packages
invisible(lapply(packages, library, character.only=T,quietly=T))



setwd('c:/users/cselinger/Dropbox (IDM)/HIV-GWAS/')#add your dropbox path....

load('Sulggi/code/cd4Setupsleigen101116.rdata')
data<-as.data.table(pheno)
data<-data[cohort=="UARTO",list(id,age,gender,cd4_num,treatdays,treatdayspp90,treatdayspp365,timeint90,timeint90365,timeint365plus)]
data<-data[!id%in% unique(data[sqrt(cd4_num)>35,id]),]#keep only ids with sqrt cd4 <35
data<-na.omit(data)
data=data[,list(id,age,gender,sqrt(cd4_num),treatdays)];setnames(data,c('id','age','gender','y','x'))
data$id <- factor(data$id)
mat<-as.data.frame(data)
mat$id<-as.factor(mat$id);mat$gender<-as.factor(mat$gender)
mat$age<-cut(mat$age,c(seq(15,50,10),100))


data[id=='MBA1001',]


mod<-lmer(sqrt(cd4_num) ~ age + gender + (bs(treatdays,degree=1) | id)+ (bs(treatdayspp90,degree=1) | id)+ (bs(treatdayspp365,degree=1) | id), data=data,na.action = na.omit)
coefs <- summary(mod)$coef

summary(mod)


predict(mod)->foo


ggplot(data[1:800,], aes(x = treatdays,y=sqrt(cd4_num-preart_cd4))) + geom_point() + facet_wrap(~id) + geom_smooth(method = "lm", formula = y ~ splines::bs(x, 3))

ggplot(data[1:800,], aes(x = treatdays,y=sqrt(cd4_num-preart_cd4))) + geom_point() + facet_wrap(~id) + geom_smooth(method = "lm", formula = y ~ x)





mat=data[,list(id,sqrt(cd4_num),treatdays)];setnames(mat,c('id','y','x'))
mat$id <- factor(mat$id)
plot1 <- xyplot(y~x, group=id, data=mat, type="b")
plot1









######EXAMPLE
# function to set up knots
knot <- function(x, knot) {(x-knot)*(x>knot)}
knots <- function(x, knots) {
  out <- sapply(knots, function(k) knot(x, k))
  colnames(out) <- knots
  out
}



breaks <- seq(0,5)
mean.slopes <- seq(5,2,-.6)
true.slopes <- sapply(mean.slopes,function(x) rnorm(n=20,sd=3,mean=x))
true.intercepts <- runif(n=20,min=30,max=50)
mat <- matrix(nr=20*10,nc=3)
for (i in 1:20)
{
  xs <- runif(n=10,min=0,max=6)
  slopes <- true.slopes[i,]
  ys <- true.intercepts[i] + rnorm(n=10,mean=0,sd=0.9) +
    sapply(xs, function(x) sum(slopes[1:floor(x)]) + (x - floor(x)) *
             slopes[ceiling(x)])
  id <- rep(i,10)
  mat[(1+10*(i-1)):(10*i),] <- cbind(xs,ys,id)
}
mat <- as.data.frame(mat)
names(mat) <- c("x","y","id")
mat$knot <- knots(mat$x, 0:5); 
m1 <- lm(y~ age+gender+knot, data=mat)
matX <- data.frame(x=0:6); 
matX$knot <- knots(matX$x, 0:5); 




mat$knot<- knots(mat$x, c(0,90,365))
m1 <- lm(y~ age+gender+knot, data=mat)
matX<-expand.grid(x=c(0,90,365,ceiling(max(mat$x))),age=factor(levels(mat$age)[1:4]),gender=factor(levels(mat$gender)))
#matX<-data.frame(x=c(0,90,365,ceiling(max(mat$x))))
matX$knot <- knots(matX$x, c(0,90,365))


matX$predict <- predict(m1, newdata=matX)



# plot of data and predicted values
plot1 <- xyplot(y~x, group=id, data=mat, type="b")
plot2 <- xyplot(predict~x, data=matX, type="b", col="black", lwd=3,groups=c(age,gender))
plot1+plot2


ggplot(matX)+geom_line(aes(x=x,y=predict))+facet_grid(age~gender)
ggplot(mat)+geom_line(aes(x=x,y=y,color=id))+facet_grid(age~gender)+scale_color_discrete(guide = FALSE)

# piecewise curves with random effects
m2 <- lmer(y~ age+gender+knot + (knot|id), data=mat)

dfs[,numcols] <- scale(dfs[,numcols])
m4 <- update(m1,data=dfs)

summary(m2)
ranef(m2)
coef(summary(m2))[ , "Estimate"]


# get predicted values
ids <- unique(mat$id)
matXid <- expand.grid(id=unique(mat$id), x=c(0,90,365,ceiling(max(mat$x))),age=factor(levels(mat$age)[1:4]),gender=factor(levels(mat$gender)))
matXid$knot <- knots(matXid$x, c(0,90,365))
matXid$predict <- rowSums(as.matrix(coef(m2)$id[matXid$id,]) * cbind(1,matXid$age,matXid$gender,matXid$knot))

# plot of data and predicted values
plot1 <- xyplot(y~x|id, data=mat[mat$id%in%mat$id[1:25],], type="b", as.table=TRUE)
plot2 <- xyplot(predict~x|id, data=matXid[matXid$id%in%matXid$id[1:36],], type="b", col="black", lwd=1)

