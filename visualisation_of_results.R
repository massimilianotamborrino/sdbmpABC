
#----------------------------------------------------------
# Import the kept (marginal) ABC posterior samples
#----------------------------------------------------------

signew<-t(read.table("sig.txt",header=F))
munew<-t(read.table("mu.txt",header=F))
Cnew<-t(read.table("C.txt",header=F))
sigMat<-data.frame(signew)
sigkept<-sigMat$signew
muMat<-data.frame(munew)
mukept<-muMat$munew
CMat<-data.frame(Cnew)
Ckept<-CMat$Cnew

#------------------------------------------------------------
# Plot the marginal posterior densities
# and store the figure "Fig.pdf"
#------------------------------------------------------------

name<-paste("Fig.pdf",sep="")
pdf(name,width=10,height=6)
par(mfrow=c(2,3),mai = c(0.6, 0.6, 0.2, 0.03))

l<-1300
r<-2700
true_val<-2000
dens<-density(sigkept,from=l,to=r,n=10000)
plot(dens$x,dens$y,type="l",col="dodgerblue3",xaxt="n",xlab=expression(sigma),cex.lab=1.5,ylab="",xlim=c(l,r),ylim=c(0,0.005),lwd=1.5)
curve(dunif(x,min=l,max=r),add=TRUE,col="firebrick3",lwd=1)
abline(v=true_val,col="black",lwd=1)
axis(1,seq(from=1400,to=2600,by=300))
legend("topright",legend=c(expression(pi(sigma)),expression(paste(pi[ABC]^num)(paste(sigma,"|",y)))),
       lty=c(1,1),col=c("firebrick3","dodgerblue3"),cex=1.6,lwd=2,seg.len=1.0,bty="n")

l<-160
r<-280
true_val<-220
dens<-density(mukept,from=l,to=r,n=10000)
plot(dens$x,dens$y,type="l",col="dodgerblue3",xaxt="n",yaxt="n",xlab=expression(mu),cex.lab=1.5,ylab="",xlim=c(l,r),ylim=c(0,0.14),lwd=1.5)
curve(dunif(x,min=l,max=r),add=TRUE,col="firebrick3",lwd=1)
abline(v=true_val,col="black",lwd=1)
axis(1,seq(from=l,to=r,by=20))
axis(2,seq(0,0.14,0.03))
legend("topright",legend=c(expression(pi(mu)),expression(paste(pi[ABC]^num)(paste(mu,"|",y)))),
       lty=c(1,1),col=c("firebrick3","dodgerblue3"),cex=1.6,lwd=2,seg.len=1.0,bty="n")

l<-129
r<-141
true_val<-135
dens<-density(Ckept,from=l,to=r,n=10000)
plot(dens$x,dens$y,type="l",col="dodgerblue3",xaxt="n",xlab=expression(C),cex.lab=1.5,ylab="",xlim=c(l,r),ylim=c(0,0.85),lwd=1.5)
curve(dunif(x,min=l,max=r),add=TRUE,col="firebrick3",lwd=1)
abline(v=true_val,col="black",lwd=1)
axis(1,seq(from=l,to=r,by=2))
legend("topright",legend=c(expression(pi(C)),expression(paste(pi[ABC]^num)(paste(C,"|",y)))),
       lty=c(1,1),col=c("firebrick3","dodgerblue3"),cex=1.6,lwd=2,seg.len=1.0,bty="n")

plot(sigkept,mukept,type="p",col="blue",pch=19,cex=.8,xlab=expression(sigma),ylab=expression(mu),lwd=1.5,cex.lab=1.5)
abline(v=sig)
abline(h=mu)
plot(sigkept,Ckept,type="p",col="blue",pch=19,cex=.8,xlab=expression(sigma),ylab=expression(C),lwd=1.5,cex.lab=1.5)
abline(v=sig)
abline(h=C)
plot(mukept,Ckept,type="p",col="blue",pch=19,cex=.8,xlab=expression(mu),ylab=expression(C),lwd=1.5,cex.lab=1.5)
abline(v=mu)
abline(h=C)
dev.off()


