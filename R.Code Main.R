##R Script: Main

#Monthly Temperature in Dubeque Iowa, USA 2022
install.packages("TSA")
library(TSA)

#more periodic time series
par(mfrow=c(1,1))
data(tempdub) #load monthly avg temp in dubuque Iowa over 11 years
tempdub
plot(tempdub,type='o',pch=16)

#less periodic time series
data(hare)
hare
plot(hare)

#First Method to Analyze periodicity: Taken's Embedding

#select embedding dimension
d=2

#select delay parameter 
tau=10

#More Periodic Data Set
#total number data points 
Tee=length(tempdub)

#construct pt cloud vecs (2d): d=2 and tau=10 => x_t -> v_t=c(x[t],x[t+tau])'
vecs=c()
for (t in 1:(Tee-(d-1)*tau)){
  pt=c(tempdub[t],tempdub[t+(d-1)*tau])
  vecs=rbind(vecs,pt)
}
plot(vecs,type='p',xlab='t',ylab='v_t',main='pt cloud more periodic series')

#construct persistence diagram of filtered pt cloud
install.packages("TDA")
library(TDA)
diag1=ripsDiag(vecs, maxdimension=1, maxscale=max(dist(vecs)))
plot(diag1$diagram)

#Less Periodic data set
#total number data points
Tee2=length(hare)

#construct pt cloud
vecs2=c()
for (t in 1:(Tee2-(d-1)*tau)){
  pt=c(hare[t],hare[t+(d-1)*tau])
  vecs2=rbind(vecs2,pt)
}
plot(vecs2,type='p',xlab='t',ylab='v_t',main='pt cloud less periodic series')

#construct persistence diagram of filtered pt cloud
diag2=ripsDiag(vecs2, maxdimension=1, maxscale=max(dist(vecs2)))
plot(diag2$diagram)

###################

#Same Process with Synthetic Data

#construct data
tee=seq(0,2*pi,length.out=100) #input points t on [0,2*pi]

#low_freq series
low_freq=2 
f1=function(t){sin(low_freq*t)}

#mid_freq series
mid_freq=4
f2=function(t){sin(mid_freq*t)}

#high_freq series
high_freq=8
f3=function(t){sin(high_freq*t)}


#plot
par(mfrow=c(3,3))
plot(f1(tee),type='l',xlab='t',ylab='x(t)=sin(2t)')
plot(f2(tee),type='l',xlab='t',ylab='y(t)=sin(4t)')
plot(f3(tee),type='l',xlab='t',ylab='z(t)=sin(8t)')


#compute Taken's Embedding (x_t -> v_t=c(x[t],x[t+(d-1)*tau])' (2d in this case)
d=2 #select embedding dimension
tau=30 #select time lag
N=100-(d-1)*tau #number of pts in pt cloud

vecs_low=c()
for(t in 1:N){
  pt=c(f1(t),f1(t+(d-1)*tau))
  vecs_low=rbind(vecs_low,pt)
}

vecs_mid=c()
for(t in 1:N){
  pt=c(f2(t),f2(t+(d-1)*tau))
  vecs_mid=rbind(vecs_mid,pt)
}

vecs_high=c()
for(t in 1:N){
  pt=c(f3(t),f3(t+(d-1)*tau))
  vecs_high=rbind(vecs_high,pt)
}

#plot pt clouds
plot(vecs_low,type='p',xlab='',ylab='')
plot(vecs_mid,type='p',xlab='',ylab='')
plot(vecs_high,type='p',xlab='',ylab='')

#plot diagrams
install.packages("TDA")
library(TDA)

diag_low=ripsDiag(vecs_low, maxdimension=1, maxscale=max(dist(vecs_low)))
plot(diag_low$diagram,main='d=2 sin(2t)')

diag_mid=ripsDiag(vecs_mid, maxdimension=1, maxscale=max(dist(vecs_mid)))
plot(diag_mid$diagram,main='d=2 sin(4t)')

diag_high=ripsDiag(vecs_high, maxdimension=1, maxscale=max(dist(vecs_high)))
plot(diag_high$diagram,main='d=2 sin(8t)')



########################

#Comparing diagrams of dimensions 2,3,4, and 5 for x,y, and z series 
#using Wasserstein distance

#keep tau=30

#d=2: plot diagrams from code above

#compute W_{1,p}(diag_low,diag_high) #set p=1 (this is q)
wasserstein(diag_low$diagram, diag_high$diagram, p=1, dimension = 0)
    #[1] 2.74094

#compute W_{1,p}(diag_low,diag_mid) 
wasserstein(diag_low$diagram, diag_mid$diagram, p=1, dimension = 0)
    #[1] 2.371494

#compute W_{1,p}(diag_mid,diag_high) 
wasserstein(diag_mid$diagram, diag_high$diagram, p=1, dimension = 0)
    #[1] 2.408383

#compute W_{infty,p}(diag_low,diag_high) 
bottleneck(diag_low$diagram, diag_high$diagram, dimension = 0)
    #[1] 0.5001738

#compute W_{infty,p}(diag_low,diag_mid)
bottleneck(diag_low$diagram, diag_mid$diagram, dimension = 0)
    #[1] 0.2807928

#compute W_{infty,p}(diag_mid,diag_high)
bottleneck(diag_mid$diagram, diag_high$diagram, dimension = 0)
    #[1] 0.3785609


#d=3: plot persistence diagrams for d=3 (fix tau=30)
d=3 #select embedding dimension
tau=30 #select time lag
N=100-(d-1)*tau #number of pts in pt cloud

vecs_low=c()
for(t in 1:N){
  pt=c(f1(t),f1(t+tau),f1(t+(d-1)*tau))
  vecs_low=rbind(vecs_low,pt)
}

vecs_mid=c()
for(t in 1:N){
  pt=c(f2(t),f2(t+tau),f2(t+(d-1)*tau))
  vecs_mid=rbind(vecs_mid,pt)
}

vecs_high=c()
for(t in 1:N){
  pt=c(f3(t),f3(t+tau),f3(t+(d-1)*tau))
  vecs_high=rbind(vecs_high,pt)
}

diag_low=ripsDiag(vecs_low, maxdimension=1, maxscale=max(dist(vecs_low)))
plot(diag_low$diagram,main='d=3 sin(2t)')

diag_mid=ripsDiag(vecs_mid, maxdimension=1, maxscale=max(dist(vecs_mid)))
plot(diag_mid$diagram,main='d=3 sin(4t)')

diag_high=ripsDiag(vecs_high, maxdimension=1, maxscale=max(dist(vecs_high)))
plot(diag_high$diagram,main='d=3 sin(8t)')

wasserstein(diag_low$diagram, diag_high$diagram, p=1, dimension = 0)
wasserstein(diag_low$diagram, diag_mid$diagram, p=1, dimension = 0)
wasserstein(diag_mid$diagram, diag_high$diagram, p=1, dimension = 0)
bottleneck(diag_low$diagram, diag_high$diagram, dimension = 0)
bottleneck(diag_low$diagram, diag_mid$diagram, dimension = 0)
bottleneck(diag_mid$diagram, diag_high$diagram, dimension = 0)

#d=4: plot persistence diagrams for d=4 (fix tau=30)
d=4 #select embedding dimension
tau=30 #select time lag
N=100-(d-1)*tau #number of pts in pt cloud

vecs_low=c()
for(t in 1:N){
  pt=c(f1(t),f1(t+tau),f1(t+2*tau),f1(t+(d-1)*tau))
  vecs_low=rbind(vecs_low,pt)
}

vecs_mid=c()
for(t in 1:N){
  pt=c(f2(t),f2(t+tau),f2(t+2*tau),f2(t+(d-1)*tau))
  vecs_mid=rbind(vecs_mid,pt)
}

vecs_high=c()
for(t in 1:N){
  pt=c(f3(t),f3(t+tau),f3(t+2*tau),f3(t+(d-1)*tau))
  vecs_high=rbind(vecs_high,pt)
}

diag_low=ripsDiag(vecs_low, maxdimension=1, maxscale=max(dist(vecs_low)))
plot(diag_low$diagram,main='d=4 sin(2t)')

diag_mid=ripsDiag(vecs_mid, maxdimension=1, maxscale=max(dist(vecs_mid)))
plot(diag_mid$diagram,main='d=4 sin(4t)')

diag_high=ripsDiag(vecs_high, maxdimension=1, maxscale=max(dist(vecs_high)))
plot(diag_high$diagram,main='d=4 sin(8t)')

wasserstein(diag_low$diagram, diag_high$diagram, p=1, dimension = 0)
wasserstein(diag_low$diagram, diag_mid$diagram, p=1, dimension = 0)
wasserstein(diag_mid$diagram, diag_high$diagram, p=1, dimension = 0)
bottleneck(diag_low$diagram, diag_high$diagram, dimension = 0)
bottleneck(diag_low$diagram, diag_mid$diagram, dimension = 0)
bottleneck(diag_mid$diagram, diag_high$diagram, dimension = 0)







#########################

#SW1PerS Embedding & Periodicity Scoring

#recommended embedding dim and num pts in pt cloud (perea)
d=15
N=201

#re-construct original series 
tee=seq(0,2*pi,length.out=N+d) #input points t on [0,2*pi]

#low_freq series
low_freq=2 
f1=function(t){sin(low_freq*t)}

#mid_freq series
mid_freq=4
f2=function(t){sin(mid_freq*t)}

#high_freq series
high_freq=8
f3=function(t){sin(high_freq*t)}

x_t=f1(tee) #definex(t)
y_t=f2(tee) #definey(t)
z_t=f3(tee)#definez(t)


par=(mfrow=c(3,3))

#generate noise
noise=rnorm(N+d,0,0.5)
x_t.noise=f1(noise)
y_t.noise=f2(noise)
z_t.noise=f3(noise)

x_t.new=x_t+x_t.noise
y_t.new=y_t+y_t.noise
z_t.new=z_t+z_t.noise

#plot noisey series
plot(x_t.new,type='l',xlab='t',ylab='x(t)=sin(2t)',main='noisy')
plot(y_t.new,type='l',xlab='t',ylab='y(t)=sin(4t)',main='noisy')
plot(z_t.new,type='l',xlab='t',ylab='z(t)=sin(8t)',main='noisy')

#de-noise series: simple moving average window=d/3
x_t.denoised=pracma::movavg(x_t.new, d/3, type = "s") 
y_t.denoised=pracma::movavg(y_t.new, d/3, type = "s") 
z_t.denoised=pracma::movavg(z_t.new, d/3, type = "s") 

#plot denoisd series,
plot(x_t.denoised,type='l',xlab='t',ylab='x(t)=sin(2t)',main='denoised')
plot(y_t.denoised,type='l',xlab='t',ylab='y(t)=sin(4t)',main='denoised')
plot(z_t.denoised,type='l',xlab='t',ylab='z(t)=sin(8t)',main='denoised')

#Computed number of points in [0,2*pi] to evaluate g:[0,2*pi]->R
T1=N+d
  
#x(t)

#construct g: cubic spline
install.packages("purrr")
library(purrr)
g_x.t=stats::spline(seq(1,2*pi,length.out=T1), x_t.denoised, n=T1)$y

#embed time series into N d-dimensional vecs and store in x_t.vecs
x_t.vecs=plyr::ldply(map(1:N, ~g_x.t[.x:(.x+d-1)]))

#normalize/center vectos
x_t.vecs.norm=t(apply(x_t.vecs,1,FUN=function(x){(x-mean(x))/sqrt(sum((x-mean(x))^2))}))

#compute persistence diagram of normalized/centered pt.cld
diag_x.t=ripsDiag(X=x_t.vecs.norm, maxdimension = 1, maxscale = sqrt(3))

#compute score
barcodes=diag_x.t$diagram
candidates=barcodes[barcodes[,"dimension"]==1,]
max_pers_index=which.max(candidates[, 3] - candidates[, 2])
max_pers_birth=candidates[max_pers_index,2]
max_pers_death=candidates[max_pers_index,3]
S=1-((max_pers_death^2)-(max_pers_birth^2))/3

#plot diagram
plot(diag_x.t$diagram,main='S=0.822')


#y(t)

#construct g: cubic spline
install.packages("purrr")
library(purrr)
g_y.t=stats::spline(seq(1,2*pi,length.out=T1), y_t.denoised, n=T1)$y

#embed time series into N d-dimensional vecs and store in x_t.vecs
y_t.vecs=plyr::ldply(map(1:N, ~g_y.t[.x:(.x+d-1)]))

#normalize/center vectos
y_t.vecs.norm=t(apply(y_t.vecs,1,FUN=function(x){(x-mean(x))/sqrt(sum((x-mean(x))^2))}))

#compute persistence diagram of normalized/centered pt.cld
diag_y.t=ripsDiag(X=y_t.vecs.norm, maxdimension = 1, maxscale = sqrt(3))

#compute score
barcodes=diag_y.t$diagram
candidates=barcodes[barcodes[,"dimension"]==1,]
max_pers_index=which.max(candidates[, 3] - candidates[, 2])
max_pers_birth=candidates[max_pers_index,2]
max_pers_death=candidates[max_pers_index,3]
S=1-((max_pers_death^2)-(max_pers_birth^2))/3

#plot diagram
plot(diag_y.t$diagram,main='S=0.793')


#z(t)

#construct g: cubic spline
install.packages("purrr")
library(purrr)
g_z.t=stats::spline(seq(1,2*pi,length.out=T1), z_t.denoised, n=T1)$y

#embed time series into N d-dimensional vecs and store in x_t.vecs
z_t.vecs=plyr::ldply(map(1:N, ~g_z.t[.x:(.x+d-1)]))

#normalize/center vectos
z_t.vecs.norm=t(apply(z_t.vecs,1,FUN=function(x){(x-mean(x))/sqrt(sum((x-mean(x))^2))}))

#compute persistence diagram of normalized/centered pt.cld
diag_z.t=ripsDiag(X=z_t.vecs.norm, maxdimension = 1, maxscale = sqrt(3))

#compute score
barcodes=diag_z.t$diagram
candidates=barcodes[barcodes[,"dimension"]==1,]
max_pers_index=which.max(candidates[, 3] - candidates[, 2])
max_pers_birth=candidates[max_pers_index,2]
max_pers_death=candidates[max_pers_index,3]
S=1-((max_pers_death^2)-(max_pers_birth^2))/3

#plot diagram
plot(diag_x.t$diagram,main='S=0.335')



#######################

#Spectral (Fourier) Analysis of Periodicity of Time Series 
par(mfrow=c(1,1))
#construct random time series
Tee=100
t=seq(0,Tee-1,length.out=Tee)
true_freq=4/Tee #4 cycles every 40 time units
sine_vals=function(x){sin(2*pi*true_freq*x)}
x_t=sine_vals(t)
#add noise
noise=rnorm(Tee,0,0.3)
x_t.noisy=sine_vals(t)+noise

#construct spectral representation 

#num sinusoidals
q=3

#frequencies
lambda.1=1/Tee; lambda.2=1/Tee; lambda.3=4/Tee;#lambda.4=4/Tee
lambda=c(lambda.1,lambda.2,lambda.3)#,lambda.4)

#mutually-uncorrelated coefficients
var=1
a=rnorm(q,0,var)
b=rnorm(q,0,var)

#construct each term in the sum
spec.term=function(j,tee){ #j is term number
  a.j=a[j]
  b.j=b[j]
  lambda.j=lambda[j]
  return(a.j*cos(2*pi*lambda.j*tee)+b.j*sin(2*pi*lambda.j*tee))
}

#construct approximate values of x(t) using spectral representation
spec.rep.points.1=spec.term(1,t)+spec.term(2,t)+spec.term(3,t)#+spec.term(4,t)
spec.rep.points.2=spec.term(1,t)+spec.term(2,t)+spec.term(3,t)#+spec.term(4,t)
spec.rep.points.3=spec.term(1,t)+spec.term(2,t)+spec.term(3,t)#+spec.term(4,t)
spec.rep.points.4=spec.term(1,t)+spec.term(2,t)+spec.term(3,t)#+spec.term(4,t)
min.spec=min(spec.rep.points.1,spec.rep.points.2,spec.rep.points.3)#,spec.rep.points.4)
max.spec=max(spec.rep.points.1,spec.rep.points.2,spec.rep.points.3)#,spec.rep.points.4)

#plot
plot(x_t.noisy,type='p',pch=16,xlab='t',ylab='x(t)',ylim=c(min.spec,max.spec),main='time series with 4 different spectral approximations')
lines(x_t.noisy,type='l',lwd=2)
lines(spec.rep.points.1,col='blue',lwd=2)
lines(spec.rep.points.2,col='green',lwd=2)
lines(spec.rep.points.3,col='orange',lwd=2)
lines(spec.rep.points.4,col='yellow',lwd=2)
legend("bottomright", legend=c("x.t", "spec.rep.1", "spec.rep.2", "spec.rep.3", "spec.rep.4"), 
       col=c("black", "blue", "green", "orange", "yellow"), lty=c(1, 1, 1, 1, 1), lwd=c(2, 2, 2, 2, 2), pch=c(16, NA, NA, NA, NA))

###################
#plot a periodogram of series against various values of frequency 
j=seq(1,Tee/2,length.out=Tee/2)
lambda=j/Tee

#cosine transform
C.lambda=c()
for (i in 1:length(j)){
  C.lambda.i=(1/sqrt(Tee))*sum(x_t.noisy*cos(2*pi*lambda[i]*t))
  C.lambda[i]=C.lambda.i
}

#sine transform
S.lambda=c()
for (i in 1:length(j)){
  S.lambda.i=(1/sqrt(Tee))*sum(x_t.noisy*sin(2*pi*lambda[i]*t))
  S.lambda[i]=S.lambda.i
}

#periodogram vals
I.lambda=C.lambda^2+S.lambda^2

#plot periodogram against frequency-val
plot(lambda,I.lambda,xlab='lambda',ylab='I(lambda)',main='periodogram',type='l')

#find lambda coresponding to max I(lambda.) value..
best.lambda=lambda[which.max(I.lambda)]
best.lambda

#best lambda is 0.04 which is true frequency
abline(v=best.lambda,col='blue')
legend("topright",legend='frequency approx = 0.04',col='blue',lty=1,lwd=1,pch=NA)

#compute persistence diagram of periodogram (2nd ord spectra)
diag.x=gridDiag(FUNvalues=I.lambda,location = FALSE, sublevel=TRUE)$diagram
plot(diag.x)
title(main=c('Diagram from Sublevel Set', 'Filtration on Periodogram'))

######################

#Compute Pers Diag of DTM of a time series 
library(TDA)
par(mfrow=c(2,2))#, mar=c(4,4,2,2)+.1, oma=c(2,2,0,0))

#use z(t)=sin(2*pi*8*t) (frequency 8 for every 2*pi units)

#construct & plot the time series
Tee=400
t=seq(0,Tee-1,length.out=Tee)
z.t=function(x){
  sin(2*pi*(8/Tee)*x)
}
series=z.t(t)
plot(t,series,type='l',xlab='',ylab='')

#construct and plot the point cloud 

#Taken's embedding 
d=2
tau=5
N=Tee-(d-1)*tau
PC=c()
for (i in 1:N){
  PC=rbind(v.t,c(series[i],series[i+(d-1)*tau]))
}
plot(PC,type='p',xlab='',ylab='')

#compute DTM of z(t) 

#define tuning param m in (0,1)
m0=0.05; 

#define step size for values of t in Grid
by=0.05; 

#define range of time (t) and sequence (z(t)) values in point cloud as xlim and ylim, respectively
Xlim=range(PC[,1])
Ylim=range(PC[,2])

#construct (t,z(t))-grid pts to plug into DTM function
Xseq=seq(from=Xlim[1],to=Xlim[2],by=by)
Yseq=seq(from=Ylim[1],to=Ylim[2],by=by)
Grid=expand.grid(Xseq,Yseq);

#compute DTM function on grid points with tuning param m=0.05
DTM=dtm(X=PC,Grid=Grid,m0=m0)

#plot DTM points
  #z is a matrix of DTM outputs for every grid pt (t,z(t)) (t=row, z(t)=column)
persp(x = Xseq, y = Yseq, z = matrix(DTM, nrow = length(Xseq), ncol = length(Yseq)),
      xlab = "t", ylab = "z(t)", zlab = "DTM(t,z(t))", theta = -20, phi = 35, scale = FALSE,
      expand = 4, col = "red", border = NA, ltheta = 10, shade = 0.5,asp=1)

#plot persistence diagram resulting from sublevel set filtration on the DTM function for t in [0,infty)
DTM.diag=gridDiag(X=PC, FUN=dtm, lim=cbind(Xlim, Ylim), by=by, m0=m0) 
plot(DTM.diag$diagram)

#visualize sublevel set L_{0.30}={(x,y) in DTM.vecs;z<=0.30} of DTM function
par(mfrow=c(1,3))
#we want to plot all (y,z) such that x=0.3
z = matrix(DTM, nrow = length(Xseq)); z=as.vector(z)
DTM.vecs=cbind(Grid[,1],Grid[,2],z)
L_0.1=DTM.vecs[DTM.vecs[,3]<=0.1,]
plot(L_0.1[,1:2],type='p',pch=16,col='red',main='L_{0.10}(DTM)',xlab='y',ylab='z')
L_0.15=DTM.vecs[DTM.vecs[,3]<=0.15,]
plot(L_0.15[,1:2],type='p',pch=16,col='red',main='L_{0.15}(DTM)',xlab='y',ylab='z')
L_0.5=DTM.vecs[DTM.vecs[,3]<=0.5,]
plot(L_0.5[,1:2],type='p',pch=16,col='red',main='L_{0.5}(DTM)',xlab='y',ylab='z')


###################

#Walsh Transforms & Plotting Periodograms 
par(mfrow=c(3,2))

Tee=400
t=seq(0,Tee-1,length.out=Tee)
true_frequency=4/400 #num cycles divided by total num pts
x=function(tee){sin(2*pi*true_frequency*tee)}
x.t=x(t)
noise=rnorm(400,0,0.3)
x.t_noisy=x.t+noise
plot(t,x.t_noisy,type='l',ylab='x(t)',main='')

walsh.f=function(n,Tee){
  lambda=seq(0,1,length.out=Tee)
  walsh.vals=c()
  for (k in 0:(n-1)){
    switch.interval.length=1/n
    interval.min=k*switch.interval.length
    interval.max=(k+1)*switch.interval.length
    walsh.indices=which(lambda>=interval.min & lambda<=interval.max)
    walsh.vals[walsh.indices]=(-1)^k
  }
  return(walsh.vals)
}

#generate/plott W(3,lambda) for Tee=400 pts lambda in [0,1]
Walsh.3=walsh.f(8,400)
lines(t,Walsh.3,col='blue',lwd=2)
legend("topright",legend=c('x(t)','W(7,lambda)'),col=c('black','blue'),lty=1,lwd=1,pch=NA,cex=0.7)

walsh.f.mod=function(n,lambda){
  walsh.vals=c()
  for (k in 0:(n-1)){
    switch.interval.length=1/n
    interval.min=k*switch.interval.length
    interval.max=(k+1)*switch.interval.length
    walsh.indices=which(lambda>=interval.min & lambda<=interval.max)
    walsh.vals[walsh.indices]=(-1)^k
  }
  return(walsh.vals)
}

#generate a Walsh periodogram
Tee=400
j=seq(1,Tee-1,length.out=Tee-1)
lambda=j/Tee
periodogram.pts=c()
for (i in 1:length(lambda)){
  lambda.j=lambda[i]
  Walsh.t=c()
  for (t in 1:Tee){
    Walsh.t[t]=walsh.f.mod(t,lambda.j)
  }
  periodogram.pts[i]=((1/sqrt(Tee))*sum(x.t*Walsh.t))^2
}

plot(lambda,periodogram.pts,type='l',xlab='switch rate',ylab='',main='Walsh Periodogram',xlim=c(0,0.1))
best.lambda=lambda[which.max(periodogram.pts)]
abline(v=best.lambda,col='blue')
legend("topright",legend='sequency approx = 0.02',col='blue',lty=1,lwd=1,pch=NA,cex=0.7)

walsh.diag=gridDiag(FUNvalues=periodogram.pts,location = FALSE, sublevel=TRUE)
plot(walsh.diag$diagram,main='Diagram After Sublevel Set Filtration')

###################

#computing persistence diagrams of modified walsh-fourier transforms of time series
#using sublevel set filtration on smooth tapered walsh-fourier periodogram

#plot original time series
par(mfrow=c(3,3))
Tee=400
t=seq(0,Tee-1,length.out=Tee)
x=function(tee){sin(2*pi*(2/400)*tee)}
x.t=x(t)
plot(t,x.t,type='l',main='frequency=2/400')
y=function(tee){sin(2*pi*(4/400)*tee)}
y.t=y(t)
plot(t,y.t,type='l',main='frequency=4/400')
z=function(tee){sin(2*pi*(8/400)*tee)}
z.t=z(t)
plot(t,z.t,type='l',main='frequency=8/400')

#compute smooth tapered periodograms using regression model x(t)=b0+b1*t+e(t) 
residuals=lm(x.t~t)$residuals #linear model for series is x(t)=0.48-0.002*t+e(t)
sigma=sd(residuals)
ts.x=residuals/sigma
smooth.periodogram.x=spec.pgram(ts.x, kernel=kernel("modified.daniell", c(1)), plot = F)
periodogram.x.pts=smooth.periodogram.x$spec
plot(seq(1,length(periodogram.x.pts),length.out=length(periodogram.x.pts)),periodogram.x.pts,type='l',xlab='',ylab='smooth periodogram',main='switch approx=2',xlim=c(0,10))

#compute smooth tapered periodograms using regression model y(t)=b0+b1*t+e(t) 
residuals=lm(y.t~t)$residuals 
sigma=sd(residuals)
ts.y=residuals/sigma
smooth.periodogram.y=spec.pgram(ts.y, kernel=kernel("modified.daniell", c(1)), plot = F)
periodogram.y.pts=smooth.periodogram.y$spec
plot(seq(1,length(periodogram.y.pts),length.out=length(periodogram.y.pts)),periodogram.y.pts,type='l',xlab='',ylab='smooth periodogram',main='switch approx=4',xlim=c(0,10))

#compute smooth tapered periodograms using regression model z(t)=b0+b1*t+e(t) 
residuals=lm(z.t~t)$residuals 
sigma=sd(residuals)
ts.z=residuals/sigma
smooth.periodogram.z=spec.pgram(ts.z, kernel=kernel("modified.daniell", c(1)), plot = F)
periodogram.z.pts=smooth.periodogram.z$spec
plot(seq(1,length(periodogram.z.pts),length.out=length(periodogram.z.pts)),periodogram.z.pts,type='l',xlab='',ylab='smooth periodogram',main='switch approx=8',xlim=c(0,10))


#compute persistence diagrams of smooth tapered periodograms using sublevel set filtration
diag.tsx=gridDiag(FUNvalues=smooth.periodogram.x$spec,location = FALSE, sublevel=TRUE)$diagram
diag.tsz=gridDiag(FUNvalues=smooth.periodogram.z$spec,location = FALSE, sublevel=TRUE)$diagram
diag.tsy=gridDiag(FUNvalues=smooth.periodogram.y$spec,location = FALSE, sublevel=TRUE)$diagram
plot(diag.tsx)
plot(diag.tsy)
plot(diag.tsz)

######### Review of Walsh Transforms


#W(t,j) for j=0,..,T-1
  #interpretation:
    #let t=0,..,T-1 and k=0,..,t
    #then the walsh function with t switches over every T time units is equal to
    #(-1)^k for indices j=k*(T/(t+1)),..,((k+1)*(T/(t+1)))-1 

#define general definition of W(t,j) for all t and j
Walsh.t.j=function(tee,Tee_val,jay){
  walsh.vals=c()
  for (k in 0:tee){
    j.start=round(k*(Tee_val/(tee+1)))
    j.end=round((k+1)*((Tee_val/(tee+1)))-1)
    walsh.vals[j.start:j.end]=(-1)^k
  }
  if (j.end!=length(jay)){walsh.vals[j.end:length(jay)]=(-1)^k}
  return(walsh.vals)
}

#plot noisy x(t) (sin curve with frequency 2 cycles every T=256 time pts) 
  ###and W(t,j) for t=0,1,2,3 using defined walsh function

#define noisy x(t)
Tee=2^8 #suppose x(t) is a sequence with 2^8=256 measures
t=seq(0,Tee-1,length.out=Tee) #define the Tee evenly-spaced time indices from 1 to Tee
lambda.j=j/Tee #define sequencies (switch-rates)
freq=2/Tee #2 cycles per every period of Tee time units
x=function(tee){sin(2*pi*(freq)*tee)} #define series x(t) with frequency 2/256
noise=rnorm(Tee,0,0.2) #define noise
x.t=x(t)+noise #define noisy series

#define sequency (switch-rate) indices
j=seq(0,Tee-1,length.out=Tee) 

#plot x(t) with W(0,j) 
W.0.j=Walsh.t.j(0,Tee,j)
plot(t,x.t,,type='l',main='W(0,j)') 
lines(W.0.j,col='blue',type='l',lwd=2)

#plot x(t) with W(1,j)
W.1.j=Walsh.t.j(1,Tee,j)
plot(t,x.t,,type='l',main='W(1,j)') 
lines(W.1.j,col='blue',type='l',lwd=2)

#plot x(t) with W(2,j)
W.2.j=Walsh.t.j(2,Tee,j)
plot(t,x.t,,type='l',main='W(2,j)') 
lines(W.2.j,col='blue',type='l',lwd=2)

#plot x(t) with W(3,j)
W.3.j=Walsh.t.j(3,Tee,j)
plot(t,x.t,,type='l',main='W(3,j)') 
lines(W.3.j,col='blue',type='l',lwd=2)

#walsh function to output values at single index j rather than at a range of values
Walsh.val=function(tee,Tee_val,j.index){
  walsh.val=0 #default walsh value 
  #compute intervals for walsh values of \pm 1
  for (k in 0:tee){
    j.start=round(k*(Tee_val/(tee+1)))
    j.end=round((k+1)*((Tee_val/(tee+1)))-1)
    #return W(tee,j.index)=(-1)^k if j.index is in one of these intervals
    if (j.index>=j.start && j.index<=j.end){
      walsh.val=(-1)^k
      break #stop looking for W(t,j)
    }
  }
  return(walsh.val)
}


#Plot Periodogram I_{w}(j) for j=0,..,T-1
  #I(j)=[sum_{t=0}^{T-1}x(t)W(t,j)]^{2} for j=0,..,T-1
period.pts=c()
for (j.index in 1:length(j)){
  sum_walsh.terms=0 
  for (t.index in 1:Tee){
    new_walsh.term=x.t[t.index]*Walsh.val(t.index,Tee,j.index)
    sum_walsh.terms=sum_walsh.terms+new_walsh.term
  }
  period.pts[j.index]=sum_walsh.terms^2
}
plot(j,period.pts,type='l',main='Walsh Periodogram',xlim=c(0,100))

#find and add best approximation j (closest # switches to true fruequency of x)
best.j.index=which.max(period.pts)
best.j=j[best.j.index]
abline(v=best.j,col='blue')
legend("topright",legend='num zero-crossings approx = 2',col='blue',lty=1,lwd=1,pch=NA,cex=0.7)


#plot persistence diagram using sublevel set filtration on periodogram
#use gridDiag
library(TDA)
diag.x=gridDiag(FUNvalues=period.pts,location = FALSE, sublevel=TRUE)$diagram
plot(diag.x)
title(main = c("Diagram from Sublevel Set", "Filtration on Periodogram"))

#exmaple of square-wave series
t=seq(0,200,length.out=200)
j=seq(0,200,length.out=200)
Tee=200
x.t=Walsh.t.j(3,Tee,j)
noise=rnorm(200,0,0.4)
x.t_noisy=x.t+noise
plot(j,x.t_noisy,type='l',xlab='t',ylab='x(t)')


#########################################

#Persistence Landscapes : Compute for x(t), y(t), z(t)
par(mfrow=c(1,2))

#x(t)
Tee=200
t=seq(1,Tee,length.out=Tee)
low_freq=(8/Tee)
z=function(tee){sin(2*pi*low_freq*tee)}
z.t=z(t)
plot(t,z.t,type='l',main='Time Series')

#taken's embedding
d=2
tau=160
N=Tee-(d-1)*tau
vecs_low=c()
for(t in 1:N){
  pt=c(z.t[t],z.t[t+(d-1)*tau])
  vecs_low=rbind(vecs_low,pt)
}
dim(vecs_low)
plot(vecs_low,type='p',xlab='',ylab='',main='Takens Embedding')

#generate diag
diag.z=ripsDiag(vecs_low, maxdimension = 1, maxscale = max(dist(vecs_low)))
#plot(diag.z$diagram,main='Persistence Diagram')

#plot only H1 points in diagram
diag.z.1=diag.z$diagram[diag.z$diagram[,1]==1,]
x=c(diag.z.1[,2])
y=c(diag.z.1[,3])
plot(x[2:4],y[2:4],type='p',col='red',pch=16,xlim=c(0,1),ylim=c(0,1),xlab='birth',ylab='death',main='H1 Persistence Diagram')
points(x[1],y[1],col='blue',pch=16)
lines(seq(0,1,length.out=100),seq(0,1,length.out=100),lty=2,col='black')

#add 1st order pts
first.order.pt=c(diag.z.1[1,2],diag.z.1[1,3])
vertical_segment=segments(first.order.pt[1],first.order.pt[1],first.order.pt[1],first.order.pt[2],col='blue',lty=2)
horizontal_segment=segments(first.order.pt[1],first.order.pt[2],first.order.pt[2],first.order.pt[2],col='blue',lty=2)

#add 2nd order pts
second.order.pts=cbind(diag.z.1[2:4,2],diag.z.1[2:4,3])
vertical_segment.1=segments(second.order.pts[1,1],second.order.pts[1,1],second.order.pts[1,1],second.order.pts[1,2],col='red',lty=2)
horizontal_segment.1=segments(second.order.pts[1,1],second.order.pts[1,2],second.order.pts[1,2],second.order.pts[1,2],col='red',lty=2)
vertical_segment.2=segments(second.order.pts[2,1],second.order.pts[2,1],second.order.pts[2,1],second.order.pts[2,2],col='red',lty=2)
horizontal_segment.2=segments(second.order.pts[2,1],second.order.pts[2,2],second.order.pts[2,2],second.order.pts[2,2],col='red',lty=2)
vertical_segment.3=segments(second.order.pts[3,1],second.order.pts[3,1],second.order.pts[3,1],second.order.pts[3,2],col='red',lty=2)
horizontal_segment.3=segments(second.order.pts[3,1],second.order.pts[3,2],second.order.pts[3,2],second.order.pts[3,2],col='red',lty=2)

legend("bottomright",legend=c('1st order peak','2nd order peaks'),col=c('blue','red'),lty=c(2,2),lwd=c(1,1),pch=c(NA,NA),cex=c(0.7,0.7))

#persistence landscape computation of 1-th homology groups (x4)

#(1): extract 1st homology groups from diag.z
diag.z.1=diag.z$diagram[diag.z$diagram[,1]==1,]
num_1.groups=nrow(diag.z.1)
#select grid for l parameter 
#construct input values of l
M1=min(diag.z.1[,2])
M2=max(diag.z.1[,3])
num_gridpts=500
delta=(M2-M1)/num_gridpts
l.grid=c()
for (k in 0:num_gridpts){
  l.grid[k+1]=M1+k*delta
}

#(2): for each value of l in the grid, compute the 
#persistence landscape across all 1-th homology groups
#PL_{1,k}(l)={min(l-la_{1,k,1},la_{1,k,2}-l)|(la_{1,k,1},la_{1,k,2}) is the barcode of the kth 1th-homology group}
PL.1=matrix(0,ncol=num_1.groups,nrow=num_gridpts+1)
for (i in 1:(num_gridpts+1)){
  l=l.grid[i]
  for (j in 1:num_1.groups){
    birth=diag.z.1[j,2]
    death=diag.z.1[j,3]
    PL.1[i,j]=ifelse(l-birth>0 & death-l>0,min(l-birth,death-l),0)
  }
}

#(3): fix l (row) and sort each persistence landscape density (column entry) in decreasing order (for each l)
for (i in 1:(num_gridpts+1)){
  PL.1[i,]=sort(PL.1[i,],decreasing=TRUE)
}

#(4): for v=0,..,num_1.groups (column entry), output the vth order persistence landscape PL_{1,v}(l) for each l (row),
#where PL_{1,v}={vth entry of lth row if vth entry is >0, 0 otherwise}

#plot vth order persistence landscape of 1-th homology
#groups for v=1,2,3,4

#v=1
plot(l.grid,PL.1[,1],lty=2,type='l',pch=16,col='blue',xlab='l',ylab='PL.v(l)'
     ,xlim=c(min(l.grid),max(l.grid)),ylim=c(min(PL.1),max(PL.1)))
title(c("Vth Order Persistence Landscapes", 
         "of 1-Homology Groups"))

#v=2
lines(l.grid,PL.1[,2],lty=2,col='red')

#v=3
lines(l.grid,PL.1[,3],lty=2,col='green')

#v=4
lines(l.grid,PL.1[,4],lty=2,col='purple')

legend("topright",legend=c('v=1','v=2','v=3','v=4'),col=c('blue','red','green','purple'),lty=c(2,2,2,2),lwd=c(1,1,1,1),pch=c(NA,NA,NA,NA),cex=c(0.7,0.7,0.7,0.7))


#find diagram pts. in persistence landscapes

#v=1: Find maximum (x1) of blue PL pts.
PL.1.max.index=which.max(PL.1[,1])
PL.1.max=PL.1[,1][PL.1.max.index]
max.v1=c(l.grid[PL.1.max.index],PL.1.max)
points(max.v1[1],max.v1[2],col='blue',pch=16)

#v=2: Find maximums (x3) of red PL pts.
library(pracma)
PL.1.max.index=findpeaks(PL.1[,2])
max.v1.1=c(l.grid[PL.1.max.index[1,2]],PL.1.max.index[1,1])
max.v1.2=c(l.grid[PL.1.max.index[2,2]],PL.1.max.index[2,1])
max.v1.3=c(l.grid[PL.1.max.index[3,2]],PL.1.max.index[3,1])
points(max.v1.1[1],max.v1.1[2],col='red',pch=16)
points(max.v1.2[1],max.v1.2[2],col='red',pch=16)
points(max.v1.3[1],max.v1.3[2],col='red',pch=16)

#v=3: Find maximums (x1) of green PL pts.
PL.1.max.index=which.max(PL.1[,3])
PL.1.max=PL.1[,3][PL.1.max.index]
max.v1.3=c(l.grid[PL.1.max.index],PL.1.max)
points(max.v1.3[1],max.v1.3[2],col='green',pch=16)

#v=4: no max

#persistence landscape computation of 0-th homology groups (x4)

#(1): extract 1st homology groups from diag.z
diag.z.0=diag.z$diagram[diag.z$diagram[,1]==0,]
num_0.groups=nrow(diag.z.0)
#select grid for l parameter 
#min M1
M1=min(diag.z.0[,2])
M2=max(diag.z.0[,3])
num_gridpts=500
delta=(M2-M1)/num_gridpts
l.grid=c()
for (k in 1:num_gridpts){
  l.grid[k]=M1+k*delta
}

#(2): for each value of l in the grid, compute the pe
#persistence landscape across all 1-th homology groups
#PL_{1,k}(l)={min(l-la_{1,k,1},la_{1,k,2}-l)|(la_{1,k,1},la_{1,k,2}) is the barcode of the kth 1th-homology group}
PL.0=matrix(,ncol=num_0.groups,nrow=num_gridpts)
for (i in 1:num_gridpts){
  l=l.grid[i]
  for (j in 1:num_0.groups){
    birth=diag.z.0[j,2]
    death=diag.z.0[j,3]
    PL.0[i,j]=ifelse(l-birth>0 & death-l>0,min(l-birth,death-l),0)
  }
}

#(3): fix l (row) and sort each persistence landscape density (column entry) in decreasing order (for each l)
for (i in 1:num_gridpts){
  PL.0[i,]=sort(PL.0[i,],decreasing=TRUE)
}

#(4): for v=0,..,num_1.groups (column entry), output the vth order persistence landscape PL_{1,v}(l) for each l (row),
#where PL_{1,v}={vth entry of lth row if vth entry is >0, 0 otherwise}
v.orders=1:num_0.groups
PL.0.v=matrix(,ncol=num_0.groups,nrow=num_gridpts)
for (v in 1:length(v.orders)){
  for (l.index in 1:num_gridpts){
    PL.0.v[l.index,v]=ifelse(PL.0[l.index,v]>0,PL.0[l.index,v],0)
  }
}

#(5): plot vth order persistence landscape of 0-th homology
#groups for v=1,..,40

#v=1
plot(1:num_gridpts,PL.0.v[,1],lty=2,type='l',col='blue',xlab='Distance',ylab='Density')
title(c("0-th Homology Groups",
      "40 Total"))
#v=5
lines(1:num_gridpts,PL.0.v[,5],lty=2,col='red')

#v=10
lines(1:num_gridpts,PL.0.v[,10],lty=2,col='green')

#v=20
lines(1:num_gridpts,PL.0.v[,20],lty=2,col='purple')

#v=40
lines(1:num_gridpts,PL.0.v[,40],lty=2,col='orange')

legend("topright",legend=c('v=1','v=5','v=10','v=20','v=40'),col=c('blue','red','green','purple','orange'),lty=c(2,2,2,2,2),lwd=c(1,1,1,1,1),pch=c(NA,NA,NA,NA,NA),cex=c(0.7,0.7,0.7,0.7,0.7))



###################################


#Applications of Persistence Landscapes

#Construct Wang et. al's (2018) smoothed weighted fourier transform as a morse function
#to compute persistence landscapes of all orders
par(mfrow=c(2,2))

#we apply it to our original series y(t) with frequency 4
Tee=500
t=seq(1,Tee,length.out=Tee)
freq=4/Tee
y=function(tee){sin(2*pi*freq*tee)}
y.t=y(t)
plot(t,y.t,type='l',main='',ylim=c(min(min(PL.vals,y.t)),max(max(PL.vals,y.t))),lwd=1)

#select max number of cycles k and number of pts per period n
k=99
n=500

#construct varying cycle numbers j for fourier sum
j=0:k

#construct weights for cosine in fourier sum a_{j} for j=0,..,k
a.j=c()
for (j.index in 1:k){
  fourier.term=c()
  for (t.index in 1:Tee){
    fourier.term[t.index]=y.t[t.index]*cos(((2*pi*j[j.index])/Tee)*t[t.index])
  }
  a.j[j.index]=(2/Tee)*sum(fourier.term)
}
#add a_{0}
a.0=mean(y.t); a.j=c(a.0,a.j) 

#construct weights for sine in fourier sum b_{j} for j=1,..,k
b.j=c()
for (j.index in 1:k){
  fourier.term=c()
  for (t.index in 1:Tee){
    fourier.term[t.index]=y.t[t.index]*sin(((2*pi*j[j.index])/Tee)*t[t.index])
  }
  b.j[j.index]=(2/Tee)*sum(fourier.term)
}

#compute standard deviation (user chosen) 
sigma=1

#compute s
a.m=median(abs(a.j[2:k+1]))
b.m=median(abs(b.j))
s=median(c(abs(a.j[2:k+1]-a.m),abs(b.j-b.m)))

#compute T_{u}
T.u=s*sqrt(2*log(n))

#compute I_{1} and I_{2}
I.1=which(abs(a.j)>T.u)
I.2=which(abs(b.j)>T.u)

#compute smoothed weighted fourier transform

#construct morse function
mu.Tu.k=function(tee){sum(exp(-(2*pi*I.1/Tee)^2)*sigma*a.j[I.1]*cos((2*pi*I.1/Tee)*tee)) +
  sum(exp(-(2*pi*I.2/Tee)^2)*sigma*b.j[I.2]*sin((2*pi*I.2/Tee)*tee))}

#output transform values mu_{T.u}^{k}(t) for all t=1,..,Tee
mu.Tu.k.vals=c()
for (t.index in 1:Tee){
  mu.Tu.k.vals[t.index]=mu.Tu.k(t[t.index])
}

#plot smoothed weighted fourier transform (morse function)
lines(t,mu.Tu.k.vals,col='blue')
title(c("Smoothed Fourier
      Transform of y(t)"))
legend("topright",legend=c('y(t)','morse(t)'),col=c('black','blue'),lty=c(1,1))

#persistence diagram of sublevel set filtration of morse function
#diag.mu=gridDiag(FUNvalues=mu.Tu.k.vals,location = FALSE, sublevel=TRUE)$diagram
#plot(diag.mu)

#persistence diagram of Taken's Embedding of morse function (Wang et. al used sublevel set filtration instead)
d=2
tau=23
N=Tee-(d-1)*tau
vecs=c()
for(t in 1:N){
  pt=c(mu.Tu.k.vals[t],mu.Tu.k.vals[t+(d-1)*tau])
  vecs=rbind(vecs,pt)
}
dim(vecs)
plot(vecs,type='p',xlab='',ylab='',pch=16)
title(c('Takens Embedding of
       Transformed Series'))
diag.mu=ripsDiag(vecs, maxdimension = 1, maxscale = max(dist(vecs)))$diagram
plot(diag.mu)
title(c('Diagram After VR
        Filtration on Embedded Pts'))
                                                         
#plot v-th order landscapes of 1-homology groups of diagram 
#(Wang et al used sublevel set filtration but we are using VR filtration)
#Use code landscape() function as alternative: (borrowed from Ravishankar & Chen (2019))

#Instead we use our function
#(1): extract 1st homology groups from diag.z
diag.mu.1=diag.mu[diag.mu[,1]==1,]
num_1.groups=nrow(diag.mu.1)
M1=min(diag.mu.1[,2])
M2=max(diag.mu.1[,3])
num_gridpts=500
delta=(M2-M1)/num_gridpts
l.grid=c()
for (k in 1:num_gridpts){
  l.grid[k]=M1+k*delta
}

#(2): for each value of l in the grid, compute the 
#persistence landscape across all 1-th homology groups
PL.1=matrix(,ncol=num_1.groups,nrow=num_gridpts)
for (i in 1:num_gridpts){
  l=l.grid[i]
  for (j in 1:num_1.groups){
    birth=diag.z.0[j,2]
    death=diag.z.0[j,3]
    PL.1[i,j]=ifelse(l-birth>0 & death-l>0,min(l-birth,death-l),0)
  }
}

#(3): fix l (row) and sort each persistence landscape density (column entry) in decreasing order (for each l)
for (i in 1:num_gridpts){
  PL.1[i,]=sort(PL.1[i,],decreasing=TRUE)
}

#(4): for v=0,..,num_1.groups (column entry), output the vth order persistence landscape PL_{1,v}(l) for each l (row)
v.orders=1:num_1.groups
PL.1.v=matrix(,ncol=num_1.groups,nrow=num_gridpts)
for (v in 1:length(v.orders)){
  for (l.index in 1:num_gridpts){
    PL.1.v[l.index,v]=ifelse(PL.1[l.index,v]>0,PL.1[l.index,v],0)
  }
}

#(5): plot vth order persistence landscape of 0-th homology
#groups for v=1,..,40

#v=1
plot(1:num_gridpts,PL.1.v[,1],lty=2,type='l',col='black',xlab='l',ylab='PL(l)')
title(c('vth order Persistence Landscapes
      of 1-Homology Groups'))

#v=2
lines(1:num_gridpts,PL.1.v[,2],lty=2,col='black')

#v=3
lines(1:num_gridpts,PL.1.v[,3],lty=2,col='black')

#v=4
lines(1:num_gridpts,PL.1.v[,4],lty=2,col='black')

#v=5
lines(1:num_gridpts,PL.1.v[,5],lty=2,col='black')

#v=6
lines(1:num_gridpts,PL.1.v[,6],lty=2,col='black')

#v=7
lines(1:num_gridpts,PL.1.v[,7],lty=2,col='black')

#v=8
lines(1:num_gridpts,PL.1.v[,8],lty=2,col='black')

#v=9
lines(1:num_gridpts,PL.1.v[,9],lty=2,col='black')

#v=10
lines(1:num_gridpts,PL.1.v[,10],lty=2,col='black')

#v=11
lines(1:num_gridpts,PL.1.v[,11],lty=2,col='black')

#v=12
lines(1:num_gridpts,PL.1.v[,12],lty=2,col='black')

#v=13
lines(1:num_gridpts,PL.1.v[,13],lty=2,col='black')


#Construct Chen et al's (2019) Walsh-Fourier (non-Morse) transform of 
#categorical time series and compute PL's:

## Shanks FFT -> PL of x.1,x.2,x.3

Tee2=2^8 #Construct T2 to be a power of 2
t=1:Tee2
j=0:(Tee2-1) #construct sequence of # of zero-crossings j
lambda.j=j/Tee2 #construct sequence of j frequencies (j crossings per every Tee points)
N=3 #we select N=3 series

#define general definition of W(t,j) for all t and j (already defined)
Walsh.t.j=function(tee,Tee_val,jay){
  walsh.vals=c()
  for (k in 0:tee){
    j.start=round(k*(Tee_val/(tee+1)))
    j.end=round((k+1)*((Tee_val/(tee+1)))-1)
    walsh.vals[j.start:j.end]=(-1)^k
  }
  if (j.end!=length(jay)){walsh.vals[j.end:length(jay)]=(-1)^k}
  return(walsh.vals)
}

#walsh function to output values at single index j rather than at a range of values
Walsh.val=function(tee,Tee_val,j.index){
  walsh.val=0 #default walsh value 
  #compute intervals for walsh values of \pm 1
  for (k in 0:tee){
    j.start=round(k*(Tee_val/(tee+1)))
    j.end=round((k+1)*((Tee_val/(tee+1)))-1)
    #return W(tee,j.index)=(-1)^k if j.index is in one of these intervals
    if (j.index>=j.start && j.index<=j.end){
      walsh.val=(-1)^k
      break #stop looking for W(t,j)
    }
  }
  return(walsh.val)
}

#plot three categorical time series
par(mfrow=c(3,3))

#similar to x(t) : frequency 2 = 2 crossings
x.1=c()
x.1[1:83]=1
x.1[84:167]=-1
x.1[168:Tee2]=1
plot(1:Tee2,x.1,xlab='t',ylab='x(1,t)',type='l',main='2 switches')

#similar to y(t): freq 4=4 crossings
x.2=c()
x.2[1:26]=1
x.2[27:52]=-1
x.2[53:74]=1
x.2[75:91]=-1
x.2[92:Tee2]=1
plot(1:Tee2,x.2,xlab='t',ylab='x(2,t)',type='l',main='4 switches')


#similar to z(t) : frequency 8 = 8 crossings
x.3=c()
x.3[1:152]=1
x.3[153:162]=-1
x.3[163:176]=1
x.3[177:189]=-1
x.3[190:210]=1
x.3[211:223]=-1
x.3[224:236]=1
x.3[237:248]=-1
x.3[248:Tee2]=1
plot(1:Tee2,x.3,xlab='t',ylab='x(3,t)',type='l',main='8 switches')

#Walsh transformation of these series
W.1.la.j=c()
for (j.index in 1:Tee2){
  jay=j[j.index]
  Walsh.term=c()
  for (t.index in 1:Tee2){
    tee=t[t.index]
    Walsh.term[t.index]=x.1[t.index]*Walsh.val(tee,Tee2,j.index)
  }
  W.1.la.j[j.index]=(1/sqrt(Tee2))*sum(Walsh.term)
}
plot(1:Tee2,W.1.la.j,type='l',xlab='j',ylab='W(1,la.j)')

W.2.la.j=c()
for (j.index in 1:Tee2){
  jay=j[j.index]
  Walsh.term=c()
  for (t.index in 1:Tee2){
    tee=t[t.index]
    Walsh.term[t.index]=x.2[t.index]*Walsh.val(tee,Tee2,j.index)
  }
  W.2.la.j[j.index]=(1/sqrt(Tee2))*sum(Walsh.term)
}
plot(1:Tee2,W.2.la.j,type='l',xlab='j',ylab='W(2,la.j)')

W.3.la.j=c()
for (j.index in 1:Tee2){
  jay=j[j.index]
  Walsh.term=c()
  for (t.index in 1:Tee2){
    tee=t[t.index]
    Walsh.term[t.index]=x.3[t.index]*Walsh.val(tee,Tee2,j.index)
  }
  W.3.la.j[j.index]=(1/sqrt(Tee2))*sum(Walsh.term)
}
plot(1:Tee2,W.3.la.j,type='l',xlab='j',ylab='W(3,la.j)')


#plot PL's of Walsh Transformations

#compute W_{n,min} for n=1,..,N=3
W.1.min=min(W.1.la.j)
W.2.min=min(W.2.la.j)
W.3.min=min(W.3.la.j)

#compute W_{n,max} for n=1,..,3
W.1.max=max(W.1.la.j)
W.2.max=max(W.2.la.j)
W.3.max=max(W.3.la.j)

#compute W_{min}=min(W_{n,min}) for n=1,..,3 
W.min=min(W.1.min,W.2.min,W.3.min)

#compute W_{max}=max(W_{n,max}) for n=1,..,3
W.max=max(W.1.max,W.2.max,W.3.max)

L=500

#V_{1}(1,l)
V1=c()
for (l in 1:L){
  V1[l]=W.min+((l-1)*W.max-W.min)/(L-1) - W.1.min
}

#V_{2}(1,l)
V2=c()
for (l in 1:L){
  V2[l]=W.1.max-W.min-((l-1)*(W.max-W.min))/(L-1)
}

#PL(1,l) (first order PL of 1st series x_{1,t})
PL.1=c()
for (l in 1:L){
  PL.1[l]=min(c(V1[which(V1>0)],V2[which(V2>0)]))
}
plot(1:L,PL.1,xlab='l',ylab='PL(1,l)')

#V_{1}(2,l)
V1=c()
for (l in 1:L){
  V1[l]=W.min+((l-1)*W.max-W.min)/(L-1) - W.2.min
}

#V_{2}(2,l)
V2=c()
for (l in 1:L){
  V2[l]=W.2.max-W.min-((l-1)*(W.max-W.min))/(L-1)
}

#PL(2,l) (first order PL of 2nd series x_{2,t})
PL.1=c()
for (l in 1:L){
  PL.1[l]=min(c(V1[which(V1>0)],V2[which(V2>0)]))
}
plot(1:L,PL.1,xlab='l',ylab='PL(1,l)')

#V_{1}(3,l)
V1=c()
for (l in 1:L){
  V1[l]=W.min+((l-1)*W.max-W.min)/(L-1) - W.3.min
}

#V_{2}(3,l)
V2=c()
for (l in 1:L){
  V2[l]=W.3.max-W.min-((l-1)*(W.max-W.min))/(L-1)
}

#PL(3,l) (first order PL of 3rd series x_{3,t})
PL.1=c()
for (l in 1:L){
  PL.1[l]=min(c(V1[which(V1>0)],V2[which(V2>0)]))
}
plot(1:L,PL.1,xlab='l',ylab='PL(1,l)')

## Cooley-Tukey's FFT -> PL of x.1,x.2,x.3
#plot three categorical time series
par(mfrow=c(3,3))

#similar to x(t) : frequency 2 = 2 crossings
x.1=c()
x.1[1:167]=1
x.1[168:192]=2
x.1[193:Tee2]=0
plot(1:Tee2,x.1,xlab='t',ylab='x(1,t)',type='l',main='2 switches')

#similar to y(t): freq 4=4 crossings
x.2=c()
x.2[1:26]=1
x.2[27:52]=3
x.2[53:74]=1
x.2[75:91]=0
x.2[92:Tee2]=2
plot(1:Tee2,x.2,xlab='t',ylab='x(2,t)',type='l',main='4 switches')

#similar to z(t) : frequency 8 = 8 crossings
x.3=c()
x.3[1:152]=1
x.3[153:162]=0
x.3[163:176]=1
x.3[177:189]=2
x.3[190:210]=3
x.3[211:223]=2
x.3[224:236]=0
x.3[237:248]=0
x.3[248:Tee2]=0
plot(1:Tee2,x.3,xlab='t',ylab='x(3,t)',type='l',main='8 switches')

#apply Cooley-Tukey's FFT to x.1
W.1.transform=fft(x.1)
plot(1:Tee2,W.1.transform,type='l',xlab='j',ylab='W(1,la.j)')

#apply Cooley-Tukey's FFT to x.2
W.2.transform=fft(x.2)
plot(1:Tee2,W.2.transform,type='l',xlab='j',ylab='W(2,la.j)')

#apply Cooley-Tukey's FFT to x.3
W.3.transform=fft(x.3)
plot(1:Tee2,W.3.transform,type='l',xlab='j',ylab='W(3,la.j)')

#convert to Persistence Landscapes
  L=500; #num grid points l
  N=3; #num series 

  #compute the minimums of each transform
  W.1.min=min(Mod(W.1.transform))
  W.2.min=min(Mod(W.2.transform))
  W.3.min=min(Mod(W.3.transform))

  #compute the maximums of each transform
  W.1.max=max(Mod(W.1.transform))
  W.2.max=max(Mod(W.2.transform))
  W.3.max=max(Mod(W.3.transform))
  
  #compute the min of the minimums & the max of the maximums.
  W.min=min(W.1.min,W.2.min,W.3.min)
  W.max=max(W.1.max,W.2.max,W.3.max)
  
  #compute V1(n,l) & V2(n,l) for n=1,..,N
  
  #V1(1,l) & V2(1,l)
  V1.1=c()
  for (l in 1:L){
    V1.1[l]=W.min-(((l-1)*(W.max-W.min))/(L-1))-W.1.min
  }
  
  V2.1=c()
  for (l in 1:L){
    V2.1[l]=W.1.max-W.min-((l-1)*(W.max-W.min))/(L-1)
  }
    
  #V1(2,l)& V2(2,l)
  V1.2=c()
  for (l in 1:L){
    V1.2[l]=W.min-((l-1)*(W.max-W.min))/(L-1)-W.2.min
  }  
  
  V2.2=c()
  for (l in 1:L){
    V2.2[l]=W.2.max-W.min-((l-1)*(W.max-W.min))/(L-1)
  }
  
  #V1(3,l)& V2(3,l)
  V1.3=c()
  for (l in 1:L){
    V1.3[l]=W.min-((l-1)*(W.max-W.min))/(L-1)-W.3.min
  }  
  
  V2.3=c()
  for (l in 1:L){
    V2.3[l]=W.3.max-W.min-((l-1)*(W.max-W.min))/(L-1)
  }
  
  #compute PL(n,l) for n=1,..,N
  
  #plot PL(1,l)
  PL.1=c()
  for (l in 1:L){
   PL.1[l]=abs(min(V1.1[l],V2.1[l]))
  }
  plot(1:500,PL.1,type='l',xlab='l',ylab='PL(1,l)')

  #plot PL(2,l)
  PL.2=c()
  for (l in 1:L){
    PL.2[l]=abs(min(V1.2[l],V2.2[l]))
  }
  plot(1:500,PL.2,type='l',xlab='l',ylab='PL(2,l)')
  
  #plot PL(3,l)
  PL.3=c()
  for (l in 1:L){
    PL.3[l]=abs(min(V1.3[l],V2.3[l]))
  }
  plot(1:500,PL.3,type='l',xlab='l',ylab='PL(3,l)')
  
#stationarity and differencing
par(mfrow=c(2,3))  
  
#stationary series
variance=10
x.1=rnorm(200,0,sqrt(variance))
plot(x.1,ylim=c(-variance,variance),type='l',xlab='t',ylab='x(t)',main='stationary series')

#histograms of different window distributions
tau=36

#window 1
t.1=10
window.1=x.1[t.1:(t.1+tau)]
hist(window.1,freq=F)

#window 2
t.2=56
window.2=x.1[t.2:(t.2+tau)]
hist(window.2,freq=F)

#non-starionary series
Tee=200
noise=rnorm(Tee,0,0.5)
t=seq(0,Tee-1,length.out=Tee)
frequency=2/Tee
x.2=sin(2*pi*frequency*t)+noise
plot(x.2,type='l',xlab='t',ylab='noisy sine curve',main='non-stationary series')

#histograms of different window distributions
tau=50

#window 1
t.1=50
window.1=x.2[t.1:(t.2+tau)]
hist(window.1,freq=F)

#window 2
t.2=100
window.2=x.2[t.2:(t.2+tau)]
hist(window.2,freq=F)

#differencing of non-stationary series
par(mfrow=c(2,2))

plot(x.2,type='l',xlab='t',ylab='',main='non-stationary series')

tau=1
x.2.diff=c()
for (k in 1:(length(x.2)-1)){
  x.2.diff[k]=x.2[k+1]-x.2[k]
}

plot(x.2.diff,type='l',xlab='t',ylab='',main='lag 1 series')
acf(x.2,xlab='t', ylab='',main='ACF of non-stationary series')
acf(x.2.diff,xlab='t',ylab='',main='ACF of 1st-Differenced series')

#unit root tests to deduce order of differencing

#test if original series is stationary
install.packages("urca")
library(urca)
kpss_test=ur.kpss(x.2)
kpss_test

#p=0.8845 > 0.05 (x.2 is not stationary)
#apply 1st-order difference and run KPSS again:
kpss_test=ur.kpss(x.2.diff)
kpss_test

#p=0.0388<0.05 (x.2 is stationary so we are done)
#optimal order of differencing required is 1