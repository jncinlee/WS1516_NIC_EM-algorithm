#16'0227 Rcode for NIC Ching-Hsi Lee


#Generate 
'Figure 2'
#x_i~pois(lda_i)
#y_i~pois(lda_i*bta)
#x_1 missing
#converge iteration <5
#load Old Faithful data frame
n<-10
lda<-sample(20:1000,n);lda
bta<-0.6
bta.lda<-bta*lda
x<-rpois(n,lda)
y<-rpois(n,bta.lda)
length(x[-1])

plot(x,y)
#MLE
beta.hat<-sum(y)/sum(x);beta.hat;bta                 #should near bta
lambda.hat<-(x+y)/(beta.hat+1);lambda.hat;lda        #should near lda

#into iteration directly
lda1.ini<-1000                           #lambda initial value for lambda 1
bta.ini<-100                             #beta initial value
lda.ini<-seq(2000,n,100);lda.ini        
s<-30                                    #how many iteration
lda.itr<-c()
lda.itr<-dim(c(s,n-1))
lda1.itr<-c()
bta.itr<-c()
for (i in 1:s){
  bta.ini<-sum(y)/(lda1.ini+sum(x[-1]))
  lda1.ini<-(y[1]+lda1.ini)/(bta.ini+1)
  lda.ini<-(x[-1]+y[-1])/(bta.ini+1)        
  bta.itr[i]<-bta.ini
  lda1.itr[i]<-lda1.ini
  lda.itr[i]<-lda.ini
}

#iteration result
lda.itr
lda1.itr
bta.itr

#convergence vs iter
iteration<-seq(1,s,1)
par(mfrow=c(1,3))
plot(iteration,lda1.itr,type="l",main="Lambda1",ylab="lambda_1",lwd=2)
abline(h=lambda.hat[1],col="red",lwd=2)
plot(iteration,lda.itr,type="l",main="Lambda2",ylab='lambda_2',lwd=2)
abline(h=lambda.hat[2],col='red',lwd=2)
plot(iteration,bta.itr,type="l",main="Beta",ylab='beta',lwd=2)
abline(h=beta.hat,col="red",lwd=2)
par(mfrow=c(1,1))



#source from wikimedia
#https://commons.wikimedia.org/wiki/File:Em_old_faithful.gif

#load library for multivariate normal
library(mvtnorm)
#load Old Faithful data frame
data(faithful)

#setup grid for plotting
xpts <- seq(from=1,to=6,length.out=100)
ypts <- seq(from=40,to=100,length.out=100)

#initial parameter estimates (chosen to be deliberately bad)
theta <- list(
  tau=c(0.5,0.5),
  mu1=c(2.8,75),
  mu2=c(3.6,58),
  sigma1=matrix(c(0.8,7,7,70),ncol=2),
  sigma2=matrix(c(0.8,7,7,70),ncol=2)
)

#E step: calculates conditional probabilities for latent variables
E.step <- function(theta)
  t(apply(cbind(
    theta$tau[1] * dmvnorm(faithful,mean=theta$mu1,sigma=theta$sigma1),
    theta$tau[2] * dmvnorm(faithful,mean=theta$mu2,sigma=theta$sigma2)
  ),1,function(x) x/sum(x)))
#M step: calculates the parameter estimates which maximise Q
M.step <- function(T) list(
  tau= apply(T,2,mean),
  mu1= apply(faithful,2,weighted.mean,T[,1]),
  mu2= apply(faithful,2,weighted.mean,T[,2]),
  sigma1= cov.wt(faithful,T[,1])$cov,
  sigma2= cov.wt(faithful,T[,2])$cov)

#function to plot current data
plot.em <- function(theta){
  mixture.contour <- outer(xpts,ypts,function(x,y) {
    theta$tau[1]*dmvnorm(cbind(x,y),mean=theta$mu1,sigma=theta$sigma1) 
    + theta$tau[2]*dmvnorm(cbind(x,y),mean=theta$mu2,sigma=theta$sigma2)
  })
  contour(xpts,ypts,mixture.contour,nlevels=5,drawlabel=FALSE,col="red",xlab="Eruption time (mins)",ylab="Waiting time (mins)",main="Waiting time vs Eruption time of the Old Faithful geyser")
  points(faithful)
}

#plot initial contours
iter <- 1
png(filename=paste("em",formatC(iter,width=4,flag="0"),".png",sep=""))
plot.em(theta)
dev.off()

#Generate 
'Figure 3'
#run EM and plot
for (iter in 2:30){
  T <- E.step(theta)
  theta <- M.step(T)
  png(filename=paste("em",formatC(iter,width=4,flag="0"),".png",sep=""))
  plot.em(theta)
  dev.off()
}

#show value for each iteration
#collect theta variable
tt<-c()
dim(tt)<-c(30, 5)
for (iter in 1:30){
  T <- E.step(theta)
  theta <- M.step(T)
  tt[[iter]]<-theta
}
as.data.frame(tt)

#Generate 
'Figure 4'
seq<-seq(1,30,1)
mu1_erup<-c();mu1_wait<-c();mu2_erup<-c();mu2_wait<-c()
for(i in 1:30){
  mu1_erup[i]<-tt[[i]]$mu1[[1]]
  mu1_wait[i]<-tt[[i]]$mu1[[2]]
  mu2_erup[i]<-tt[[i]]$mu2[[1]]
  mu2_wait[i]<-tt[[i]]$mu2[[2]]
}
par(mfrow=c(1,4))
plot(seq, mu1_erup,type="l",main='First Cluster X-axis',xlab='iteration',ylab='eruption(min)',lwd=2)
abline(h=mu1_erup[30],col='blue',lwd=2)
plot(seq, mu1_wait,type="l",main='First Cluster Y-axis',xlab='iteration',ylab='waiting(min)',lwd=2)
abline(h=mu1_wait[30],col='blue',lwd=2)
plot(seq, mu2_erup,type="l",main='Second Cluster X-axis',xlab='iteration',ylab='eruption(min)',lwd=2)
abline(h=mu2_erup[30],col='red',lwd=2)
plot(seq, mu2_wait,type="l",main='Second Cluster Y-axis',xlab='iteration',ylab='waiting(min)',lwd=2)
abline(h=mu2_wait[30],col='red',lwd=2)
par(mfrow=c(1,1))



#Generate 
'Figure 567'
#source from
#https://en.wikibooks.org/wiki/Data_Mining_Algorithms_In_R/Clustering/Expectation_Maximization_%28EM%29

install.packages("mclust")
require("mclust")

#example from gyser faithful
x = faithful[,1]          # get the first column of the faithful data set
y = faithful[,2]          # get the second column of the faithful data set
plot(x,y)                 # plot the spread points before the clustering
model <- Mclust(faithful) # estimate the number of cluster (BIC), initialize (HC) and clusterize (EM)
data1 = faithful           # get the data set 
plot(model)     # plot the clustering results

#example from iris 
model1<-Mclust(data)
plot(model1)

#example for two gaussian N(1,1) N(5,1)
x1 = rnorm(n=20, mean=1, sd=1)   # get 20 normal distributed points for x axis with mean=1 and std=1 (1st class)
y1 = rnorm(n=20, mean=1, sd=1)   # get 20 normal distributed points for x axis with mean=1 and std=1 (2nd class)
x2 = rnorm(n=20, mean=7, sd=1)   # get 20 normal distributed points for x axis with mean=5 and std=1 (1st class)
y2 = rnorm(n=20, mean=7, sd=1)   # get 20 normal distributed points for x axis with mean=5 and std=1 (2nd class)
rx = range(x1,x2)                # get the axis x range
ry = range(y1,y2)                # get the axis y range
plot(x1, y1, xlim=rx, ylim=ry)   # plot the first class points
points(x2, y2)                   # plot the second class points
mix = matrix(nrow=40, ncol=2)    # create a dataframe matrix 
mix[,1] = c(x1, x2)              # insert first class points into the matrix
mix[,2] = c(y1, y2)              # insert second class points into the matrix
mixclust = Mclust(mix)           # initialize EM with hierarchical clustering, execute BIC and EM
plot(mixclust, data = mix)       # plot the two distinct clusters found

#example from Wreath
data(wreath) #data call in
wreathDefault <- mclustBIC(wreath) #default with max=9 class
wreathCustomize <- mclustBIC(wreath, G = 1:20, x = wreathDefault)
plot(wreathDefault)
plot(wreathCustomize, G = 1:20, legendArgs = list(x = "bottomleft"))
summary(wreathCustomize, wreath)
wreathBIC <- mclustBIC(wreath)
wreathBIC <- mclustBIC(wreath, G = 1:20, x = wreathBIC)
wreathModel <- summary(wreathBIC, data = wreath)
mclust2Dplot(data = wreath, what = "density", identify = TRUE, parameters = wreathModel$parameters, z = wreathModel$z)
par(mfrow=c(1,1))
plot(Mclust(wreath, G = 1:20))
?mclust