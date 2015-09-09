# For LD diameter size distribution graphs

setwd("~/R")
data<-read.csv('yourfile.csv')

par(mfrow=c(1,1))
dens_treatment1<-density(data$treatment1,bw=2,na.rm=TRUE)
dens_treatment2<-density(data$treatment2,bw=2,na.rm=TRUE)
# repeat the same for each column name 
plot(data$treatment1,col='red',xlab="x-axis description", 
     main="Graph Description", lwd=2, lty=3, xlim= c(0,40), ylim=c(0, 0.10))
lines(data$treatment2, col='blue', lwd=2, lty=2)
# repeat the same for each column name 
legend('topright',
       legend=c("treatment1", "treatment2"),
       col=c( "red", "blue"),
       pch=c(19,19), lty=c(3,3))



#To calculate the quantiles

treatment1<- quantile(data$treatment1, c( .50, .95, .90), na.rm=TRUE)
treatment2<- quantile(data$treatment2, c( .50, .95, .90), na.rm=TRUE)
# repeat the same for each column name
