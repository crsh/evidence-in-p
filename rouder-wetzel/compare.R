library('MCMCpack')

bf0 = function(t,n)
(1+t^2/(n-1))^(-n/2)

prior=function(g,n) dinvgamma(g,.5,n/2)

bf1a=function(g,t,n)
(1+g)^(-.5)*(1+t^2/((1+g)*(n-1)))^(-n/2)

joint=function(g,t,n)
prior(g,n)*bf1a(g,t,n)

bf1= function(t,n)
integrate(joint,lower=0,upper=Inf,t=t,n=n)

bf=function(t,n) bf0(t,n)/bf1(t,n)$value

bic=function(t,n) sqrt(n)*(1+t^2/(n-1))^(-n/2)

unit=function(t,n)  bf0(t,n)/bf1a(n,t,n)

#############################

target=function(t,b,n)
return((bf(t,n)-b)^2)

pdf('compare.pdf',width=12,height=6)
par(par.2pan)

v.lab=c(5,10,20,50,100,200,500,1000,2000,5000)
N=c(seq(5,20,1),seq(25,100,5),seq(200,1000,100),2000,5000)
x.lab=log(v.lab)
x=log(N)
g1=N
g2=g1
g3=g1
for (i in 1:length(N))
{
g1[i]=optimize(target,interval=c(0,10),b=1/3,n=N[i])$minimum
g2[i]=optimize(target,interval=c(0,10),b=1/10,n=N[i])$minimum
g3[i]=optimize(target,interval=c(0,10),b=1/30,n=N[i])$minimum
}
t=qt(.975,N-1)

plot(x,g1,axes=F,ylim=c(1.5,6),typ='l',xlim=log(c(5,20000)),xlab="Sample Size",ylab="Critical t-value") 
axis(1,at=x.lab,label=v.lab) 
axis(2)
lines(x,g2)
lines(x,g3)
lines(x,t)
text(log(5600),g1[length(g1)],expression(B[10]==3),adj=0)
text(log(5600),g2[length(g2)],expression(B[10]==10),adj=0)
text(log(5600),g3[length(g3)],expression(B[10]==30),adj=0)
text(log(5600),1.96,expression(alpha==.05),adj=0)
box()
mtext(side=3,'A.',cex=1.3,adj=0)



load('ej-graph.RData')
#look at dat1, one sample

one=as.matrix(dat1$one)
pair=as.matrix(dat1$paired)
two=as.matrix(dat1$two)

p=c(one[,6],pair[,8],two[,8])
effect=c(one[,4],pair[,6],two[,6])
bf=log10(c(one[,5],pair[,7],two[,7]))

good=(p>.001& p<.15)
mcol=c(rgb(1,0,0,.4),rgb(1,1,0,.4),rgb(.3,1,.8,.4),rgb(0,1,0,.4))


colCat=rep(4,sum(good))
colCat[bf[good]<log10(10)]=3
colCat[bf[good]<log10(3)]=2
colCat[bf[good]<log10(1)]=1




xtick=c((1:10)*.001,(2:10)*.01,.2)
xa=c(.001,.01,.05,.1)
ytick=c(.3,.4,.5,(1:10),20,30,40)
ya=c(.5,1,10,40)

xl=log10(c(.001,.15))
yl=log10(c(.3,40))
plot(log10(p[good]),bf[good],axes=F,typ='n',ylab=expression(paste("Bayes Factor  ",B[10])),
xlab="p-Value",xlim=xl,ylim=yl)

#rect(log10(.05),yl[1],log10(.15),yl[2],col=rgb(1,1,0,.1))
#rect(log10(.001),yl[1],log10(.01),yl[2],col=rgb(0,1,0,.2))
#rect(log10(.01),yl[1],log10(.05),yl[2],col=rgb(0,1,0,.05))
points(log10(p[good]),bf[good],pch=21,bg=mcol[1],cex=.9)
axis(1,at=log10(xtick),labels=NA)
axis(2,at=log10(ytick),labels=NA)
axis(2,at=log10(ya),label=ya)
axis(1,at=log10(xa),label=xa)
#abline(h=log10(c(1,3,10)),lty=3)
#box()
abline(v=log10(c(.05,.01)),lty=2)
mtext(side=3,'B.',adj=0,cex=1.3)
box()

dev.off()
