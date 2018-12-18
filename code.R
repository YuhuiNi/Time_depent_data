temperature=read.table("ssn_HadCET_mean.txt",header=TRUE)
size<-dim(temperature)
annualtemp<-apply(temperature[,2:5],1,mean)
t<-c(1659:2017)
annualtemp<-annualtemp[2:358]
t<-t[2:358]
lmod<-lm(annualtemp~t)
summary(lmod)$coefficients

lmod<-lm(annualtemp~t+I(t^2))
summary(lmod)$coefficients

lmod<-lm(annualtemp~t+I(t^2)+I(t^3))
summary(lmod)$coefficients

y<-resid(lmod)
library(timeDate)
library(timeSeries)
library(fBasics)
library(fUnitRoots)

adfTest(y)

par(mfrow=c(1,2))
y<-resid(lmod)
acf(y)
pacf(y)

arima(y,order=c(0,0,2))[1]
arima(y,order=c(0,0,2))[2]

plot(t,annualtemp,xlab="Time",ylab="annual temperature",type='l')
lines(t,fitted(lmod))

temp<-t(temperature[,2:5])
dim(temp)<-c(359*4,1)
temp<-temp[2:(359*4-1)]
seasontemp<-ts(temp,start=c(1659,2),frequency=4)
t<-time(seasontemp)
Q=factor(cycle(seasontemp))
reg<-lm(seasontemp~0+t+I(t^2)+I(t^3)+Q,na.action=NULL)
summary(reg)

y1=1.458*t-0.000804*(t^2)+0.0000001479*(t^3)
plot(y1,ylab="season temperature")


t<-c(1659:2017)
t<-t[2:358]
plot(t,annualtemp,type='l',xlab="Time")
lines(ksmooth(t, annualtemp, "normal", bandwidth=5),lty=1,lwd=2, col=4)
lines(ksmooth(t, annualtemp, "normal", bandwidth=30), lty=2,lwd=2, col=2)
lines(ksmooth(t, annualtemp, "normal", bandwidth=60), lty=3,lwd=2, col='green')
legend("topleft",legend=c('bandwidth=5','bandwidth=30','bandwidth=60'),col=c(4,2,'green'),lty=c(1,2,3))

x<-t
y<-annualtemp
NWSMOOTH<-function(h, y, x, z) {
    n<-length(y)
    s.hat<-rep(0, n)
    s.hat1<-rep(0, n)
    for (i in 1:n) {
        s.hat[i] <- dnorm((x[i] - z)/h)/h * y[i]
        s.hat1[i] <- dnorm((x[i] - z)/h)/h
    }
    z.hat <- sum(s.hat)/sum(s.hat1)
    return(z.hat)
}
CVRSS <- function(h, y, x) {
    cv <- NULL
    for (i in seq(x)) {
        cv[i] <- (y[i] - NWSMOOTH(h, y[-i], x[-i], x[i]))^2
    }
    mean(cv)
}
h <- seq(1, 60, by = 0.5)
cvrss.val <- rep(0, length(h))
for (i in seq(h)) {
    cvrss.val[i] <- CVRSS(h[i], y, x)
}
plot(h, cvrss.val, type = "b",xlab="bandwidth")

annualtemp=ts(annualtemp,start=1660,frequency = 1)
plot(annualtemp,type="l")
lines(lowess(annualtemp,f=0.1),col="red")

plot(annualtemp,type="l")
lines(smooth.spline(time(annualtemp), annualtemp,spar=0.7), lty=2, lwd=2, col=2)

winter<-temperature[2:358,2]
winter<-ts(winter,start=1660,frequency=1)
plot(winter,type="l",ylab="Temperature",main="Winter trend")
lines(lowess(winter,f=0.1),col="red",lty=1)
lines(smooth.spline(time(winter), winter, spar=0.7), lty=2, lwd=2, col=4)
legend("topleft",legend=c("lowess","smooth spline"),col=c(2,4),lty=c(1,2))



spring<-temperature[2:358,3]
spring<-ts(spring,start=1660,frequency=1)
plot(spring,type="l",ylab="Temperature",main="Spring trend")
lines(lowess(spring,f=0.1),col="red",lty=1)
lines(smooth.spline(time(spring), spring, spar=0.7), lty=2, lwd=2, col=4)
legend("topleft",legend=c("lowess","smooth spline"),col=c(2,4),lty=c(1,2))


summer<-temperature[2:358,4]
summer<-ts(summer,start=1660,frequency=1)
plot(summer,type="l",ylab="Temperature",main="Summer trend")
lines(lowess(summer,f=0.1),col="red",lty=1)
lines(smooth.spline(time(summer), summer, spar=0.7), lty=2, lwd=2, col=4)
legend("topleft",legend=c("lowess","smooth spline"),col=c(2,4),lty=c(1,2))


autumn<-temperature[2:358,5]
autumn<-ts(autumn,start=1660,frequency=1)
plot(autumn,type="l",ylab="Temperature",main="Autumn trend")
lines(lowess(autumn,f=0.1),col="red",lty=1)
lines(smooth.spline(time(autumn), autumn, spar=0.7), lty=2, lwd=2, col=4)
legend("topleft",legend=c("lowess","smooth spline"),col=c(2,4),lty=c(1,2))

