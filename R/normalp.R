dnormp<-function(x,mu=0,sigmap=1,p=2,log=FALSE){
if(!is.numeric(x)||!is.numeric(mu)||!is.numeric(sigmap)||!is.numeric(p)) 
stop (" Non-numeric argument to mathematical function")
if(p<1) stop("p must be at least equal to one")
if(sigmap<=0) stop("sigmap must be positive")
cost<-2*p^(1/p)*gamma(1+1/p)*sigmap
expon1<-(abs(x-mu))^p
expon2<-p*sigmap^p
dsty<-(1/cost)*exp(-expon1/expon2)
if(log==TRUE) dsty<-log(dsty)
dsty
}

pnormp<-function(q,mu=0,sigmap=1,p=2,lower.tail=TRUE,log.pr=FALSE){
if(!is.numeric(q)||!is.numeric(mu)||!is.numeric(sigmap)||!is.numeric(p)) 
stop (" Non-numeric argument to mathematical function")
if(p<1) stop("p must be at least equal to one")
if(sigmap<=0) stop("sigmap must be positive")
z<-(q-mu)/sigmap
zz<-abs(z)^p
zp<-pgamma(zz,shape=1/p,scale=p)
zp<-zp/2
zp<-ifelse(z<0,0.5-zp,0.5+zp)
if (lower.tail==FALSE) zp<-1-zp
if (log.pr==TRUE) zp<-log(zp)
zp
}

qnormp<-function(pr,mu=0,sigmap=1,p=2,lower.tail=TRUE,log.pr=FALSE){
if(!is.numeric(pr)||!is.numeric(mu)||!is.numeric(sigmap)||!is.numeric(p)) 
stop (" Non-numeric argument to mathematical function")
if(p<1) stop("p must be at least equal to one")
if(sigmap<=0) stop("sigmap must be positive")
if (log.pr==TRUE) pr<-log(pr)
if (lower.tail==FALSE) pr<-1-pr
zp<-ifelse(pr<0.5,0.5-pr,pr-0.5)
zp<-2*zp
qg<-qgamma(zp,shape=1/p,scale=p)
z<-qg^(1/p)
z<-ifelse(pr<0.5,-z,z)
q<-mu+z*sigmap
q
}


rnormp<-function(n,mu=0,sigmap=1,p=2,method=c("def","chiodi")){
if(!is.numeric(n)||!is.numeric(mu)||!is.numeric(sigmap)||!is.numeric(p)) 
stop (" Non-numeric argument to mathematical function")
if(p<1) stop("p must be at least equal to one")
if(sigmap<=0) stop("sigmap must be positive")
method <- match.arg(method)
if (method=="def"){
qg<-rgamma(n,shape=1/p,scale=p)
z<-qg^(1/p)
z<-ifelse(runif(n)<0.5,-z,z)
x<-mu+z*sigmap
}
if (method=="chiodi"){
i<-0
x<-c(rep(0,n))
while (i<n){
u<-runif(1,-1,1)
v<-runif(1,-1,1)
z<-abs(u)^p+abs(v)^(p/(p-1))
if (z<=1){
i<-i+1
x[i]<-u*(-p*log(z)/z)^(1/p)
x[i]<-mu+x[i]*sigmap
}
}
}
x
}

qqlinep<-function(y,p=2,...)
{
    y <- quantile(y[!is.na(y)],c(0.25, 0.75))
    x <- qnormp(c(0.25, 0.75),p=p)
    slope <- diff(y)/diff(x)
    int <- y[1]-slope*x[1]
    abline(int, slope, ...)
}


qqnormp<-function(y, ylim, p, main, xlab, ylab, ...) UseMethod("qqnormp")

qqnormp.default<-function(y, ylim, p=2, main="Normal of order p Q-Q plot",
            xlab="Theoretical Quantiles", ylab="Sample Quantiles",...) 
{
y<-y[!is.na(y)]
if (0 == (n<-length(y))) stop("y is empty")
if (missing(ylim)) ylim<-range(y)
x<-qnormp(ppoints(n),p=p)[order(order(y))]
plot(x,y,main=main, xlab=xlab, ylab=ylab, ylim=ylim, ...)
mtext(paste("p=",p), 3, 0.25)
invisible(list(x=x,y=y))
}


  estimatep<-function(x,mu,p=2,method=c("inverse","direct"))
  {
  if(!is.numeric(x)||!is.numeric(mu)||!is.numeric(p)) 
  stop (" Non-numeric argument to mathematical function")
  method<-match.arg(method)
  ssp<-sum(abs(x-mu)^p)/length(x)
  sp<-ssp^(1/p)
  xstand<-(x-mu)/sp
  sa<-sum(abs(xstand))
  sb<-sum(xstand*xstand)
  vi<-sqrt(length(x)*sb)/sa
  vi<-vi+((vi-1)/length(x))*5
  if (method=="inverse"){
  if (vi<1.1547005) pz<-11.5 
  else {
  zi<-(1.4142135-vi)/0.259513
  if (zi<0) pz<-1.0 
  else { 
  if (zi<0.6200052) 
  {tz<-zi/0.6200052; 
  yy<-((0.4738581*tz-0.4966873)*tz+1.0532646)*tz+1.2246159;
  pz<-1.0+tz^yy} 
  else {
  if (zi<0.7914632) 
  {tz<-(zi-0.6200052)/0.1714592;
  yy<-((0.5246979*tz-0.8167733)*tz+0.8805483)*tz+1.0859246;
  pz<-2.0+tz^yy}
  else {
  if (zi<0.8670333) 
  {tz<-(zi-0.7914632)/0.0755701;
  yy<-((0.0743092*tz-0.1269859)*tz+0.3588207)*tz+1.1227837;
  pz<-3.0+tz^yy}
  else {
  if (zi<0.9072536)
  {tz<-(zi-0.8670333)/0.0402203;
  yy<-((0.1097723*tz-0.2127039)*tz+0.3529203)*tz+1.0761256;
  pz<-4.0+tz^yy}
  else {
  if (zi<0.9314555)
  {tz<-(zi-0.9072536)/0.0242019;
  yy<-((0.0955441*tz-0.1891569)*tz+0.2961275)*tz+1.0631784;
  pz<-5.0+tz^yy}
  else {
  if (zi<0.9472072)
  {tz<-(zi-0.9314555)/0.0157518;
  yy<-((0.0862627*tz-0.1725326)*tz+0.256885)*tz+1.0540746;
  pz<-6.0+tz^yy}
  else {
  if (zi<0.9580557)
  {tz<-(zi-0.9472072)/0.0108484;
  yy<-((0.078785*tz-0.1581388)*tz+0.227011)*tz+1.04735;
  pz<-7.0+tz^yy}
  else {
  if (zi<0.9658545)
  {tz<-(zi-0.9580557)/0.0077988;
  yy<-((0.0663921*tz-0.1380841)*tz+0.2010053)*tz+1.0422984;
  pz<-8.0+tz^yy}
  else {
  if (zi<0.9716534)
  {tz<-(zi-0.9658545)/0.0057989;
  yy<-((0.0557199*tz-0.1184033)*tz+0.178176)*tz+1.038582;
  pz<-9.0+tz^yy}
  else pz<-10+(zi-.9716534)/.0283466
  }
  }
  }
  }
  }
  }
  }
  }
  }
  }
  }
  if (method=="direct"){
  fvi<-function(p) (vi-sqrt(gamma(1/p)*gamma(3/p))/gamma(2/p))*(vi-sqrt(gamma(1/p)*gamma(3/p))/gamma(2/p))
  pz<-optim(p,fvi,method="BFGS")$par
  if (pz<1.0) pz<-1.0
  if (pz>10.0) pz<-11.5
  }
  pp<-pz
  pp
  }



paramp<-function(x,p)UseMethod("paramp")

paramp.default<-function(x,p=NULL){
if(!is.numeric(x) || (!is.null(p) && !is.numeric(p))) stop (" Non-numeric argument to mathematical function")
ff<-function(Mp) sum(abs(x-Mp)^p)
  Mp<-mean(x)
if (is.null(p)){  
df<-2  
pp<-2
  iter<-0
  i<-0
  p<-estimatep(x,Mp,pp)
  while (abs(p-pp)>0.0001) {
  pp<-p
  op<-optim(Mp,ff,method="BFGS")
  Mp<-op$par
  p<-estimatep(x,Mp,pp)
  i<-i+1
  if (i==100) {iter<-1; break}
  }
  }
else {
if (p<1) stop("Value of p must be greater or equal then 1")
 df<-1  
  op<-optim(Mp,ff,method="BFGS")
  Mp<-op$par
  iter<-0
}
  if (p==1) Mp<-median(x)
  Sp<-((sum(abs(x-Mp)^p))/(length(x)-df))^(1/p)
  if (p>=11.5) {Mp<-mean(c(max(x),min(x)));   Sp<-(max(x)-min(x))/2 }
sd<-sqrt((length(x)-1)*var(x)/length(x))
mn<-mean(x)
ris<-list()
ris$mean<-mn
ris$mp<-Mp
ris$sd<-sd
ris$sp<-Sp
ris$p<-p
ris$iter<-iter
class(ris)<-"paramp"
ris
} 
 
print.paramp<-function(x,...){
dat<-c(Mean=x$mean,Mp=x$mp,Sd=x$sd,Sp=x$sp,p=x$p)
print(dat)
cat("\nno.conv =", as.logical(x$iter),"\n\n")
invisible(x)
}




kurtosis<-function(x=NULL,p,value=c("estimate","parameter")){
if (missing(p) && is.null(x)) stop("no arguments inserted")
if (!missing(p) && !is.numeric(p) || !is.null(x) && !is.numeric(x)) stop(" Non-numeric argument to mathematical function")
if (!missing(p) && p<1) stop("p must be at least 1") 
value <- match.arg(value)
 if (is.null(x)){
        vi<-sqrt(gamma(1/p)*gamma(3/p))/gamma(2/p)
        b2<-(gamma(1/p)*gamma(5/p))/(gamma(3/p))^2
        bp<-p+1
 }
 else {
 if (!is.numeric(x)) stop ("x must be a numerical vector")
    n<-length(x)
    if (value=="estimate") {
       cmp<-paramp(x)
       p<-cmp$p
       mp<-cmp$mp
     }
    else  mp<-paramp(x,p=p)$mp
    vi<-(sqrt(n*(sum((x-mp)^2))))/(sum(abs(x-mp)))
    b2<-(n*sum((x-mp)^4))/(sum((x-mp)^2))^2
  Sp<-((sum((abs(x-mp))^p))/n)^(1/p)
  S<-Sp^(2*p)
  bp<-sum((abs(x-mp))^(2*p))/(n*S)
  }
  RIS<-c(VI=vi,B2=b2,Bp=bp)
  RIS
  }

lmp<-function(formula,data,p) UseMethod("lmp")

lmp.default<-function(formula,data=list(),p=NULL){
 ll<-lm(formula=formula, data=data)
 r <- residuals(ll)
 s <- sqrt(deviance(ll)/df.residual(ll))
 hii <- lm.influence(ll)$hat
 rs <- r/(s * sqrt(1 - hii))
 cl<-match.call()
 cl[[1]]<-as.name("lmp")
 b<-vector()
 bb<-vector()
 nvar<-ll$rank
 for (i in 1:nvar) {bb[i]<-ll$coef[[i]]}
 fit<-vector()
 N<-nrow(ll$model)
 y<-ll$model[[1]]
 M<-as.matrix(ll$model)
 M[,1]<-1     
 res<-ll$residual
  f<-function(b){
  sum(abs(y-M%*%b)^pp)
  }
 if (is.null(p)){
  knp<-FALSE
  pp<-estimatep(res,mean(res),2)
  op<-optim(bb,f,method="BFGS")
  for (i in 1:nvar){b[i]<-op$par[i]}
  fit<-c(M%*%b) 
  res<-c(y-fit)
 p<-estimatep(res,mean(res),2)
 i1<-0
 iter<-0
 while (abs(pp-p)>0.0001 || abs(bb-b)>0.0001){
 pp<-p
 bb<-b
 op<-optim(bb,f,method="BFGS")
 for (i in 1:nvar){b[i]<-op$par[i]}
 fit<-c(M%*%b) 
 res<-c(y-fit)
 p<-estimatep(res,mean(res),2)
 i1<-i1+1
 if (i1==100) {iter<-1; break}
        }
        }
 else {
 knp<-TRUE
 pp<-p
 op<-optim(bb,f)
 for (i in 1:nvar){b[i]<-op$par[i]}
 fit<-c(M%*%b) 
 res<-c(y-fit)
 }
 for (i in 1:nvar){ll$coefficients[[i]]<-b[i]}
 names(res)<-as.character(1:N)
 ll$residuals<-res
 ll$call<-cl
 ll$knp<-knp
 ll$p<-p 
 names(fit)<-as.character(1:N)
 ll$fitted.values<-fit 
 if(knp==FALSE) ll$iter<-iter
 ll$rs<-rs
 class(ll)<-c("lmp","lm")
 ll
}

summary.lmp<-function(object,...){
a<-object
nvar<-a$rank
N<-nrow(a$model)
rdf<-ifelse(a$knp==FALSE,N-nvar-1,N-nvar)
r<-a$residuals
Y<-a$fitted.values
p<-a$p
iter<-a$iter
mpy<-paramp(Y,p=p)$mp
rss<-sum(abs(r)^p)
resvar<-rss/rdf
sigma<-(resvar)^(1/p)
ans <- object[c("call", "terms")]
    ans$p<-p
    ans$residuals <- r
    ans$coefficients <- a$coef
    ans$sigma <- sigma
    ans$rdf <- rdf
    ans$iter <- iter
    class(ans) <- "summary.lmp"
    ans
}
  
print.summary.lmp<-function(x,...){
a<-x
cat("\nCall:\n")
cat(paste(deparse(a$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
cat("Residuals:\n")
rq <- quantile(a$residuals)
names(rq) <- c("Min", "1Q", "Median", "3Q", "Max")
print(rq,4)
cat("\nCoefficients:\n")
print(a$coef,4)
cat("\nEstimate of p\n")
print(as.symbol(a$p))
cat("\nPower deviation of order p:",
    format(signif(a$sigma,4)))
cat("\n")
if(a$iter==1) cat("\nWarning: problems in convergence in lmp\n")
invisible(a)
}

plot.lmp<-function(x,...){
dat<-x
subc<-deparse(dat$call)
plot(dat$fitted,dat$resid,main="",xlab="Fitted values",ylab="Residuals")
abline(h=0,lty=3)
mtext("Residuals vs Fitted", 3, 0.25)
title(sub=subc)
op<-par(ask=TRUE)
qqnorm(x$rs,main="",ylab="Standardized residuals",xlab="Theoretical Quantiles")
mtext("Normal Q-Q plot", 3, 0.25)
title(sub=subc)
rs<-paramp(dat$residuals,p=dat$p)
res<-(dat$residuals-rs$mp)/rs$sp
qqnormp(res,p=rs$p,main="Normal of order p Q-Q plot",ylab="Standardized residuals",xlab="Theoretical Quantiles")
title(sub=subc)
plot(dat$fitted,(abs(res))^(1/rs$p),main="",xlab="Fitted values",ylab=expression(sqrt("Standardized residuals",p)))
mtext("Scale-Location plot", 3, 0.25)
title(sub=subc)
par(op)   
}

graphnp<-function(p=c(1,2,3,4,5),mu=0,sigmap=1,title="Normal of order p curves"){
if(!is.numeric(p)||!is.numeric(mu)||!is.numeric(sigmap))
stop (" Non-numeric argument to mathematical function")
n<-length(p)
leg<-vector()
for (i in 1:n){
if (p[i]<1) stop("There is a value of p less than 1")
}
if (n>5) stop("Please, insert only five values for p")
n1<-0
xfit<-seq(mu-4*sigmap,mu+4*sigmap,0.01)
n1<-length(p[p>=50])
if (n==1 && p>=50 || n==n1){
plot(xfit,dunif(xfit,min=mu-sigmap,max=mu+sigmap),xlim=range(mu-4*sigmap,mu+4*sigmap),
ylim=range(0,(1+0.1)/(2*sigmap)),
type="l",main=title,col=1,xlab="x",ylab="f(x)")
leg<-c(leg,expression(paste("p-> ",infinity)))
legend(mu+2.5*sigmap, (1-0.1)/(2*sigmap), leg,col=1, lty=1)
}
else {
p<-sort(p)
n<-n-n1
p<-p[1:n]
if (n1>0) { plot(xfit,dunif(xfit,min=mu-sigmap,max=mu+sigmap),xlim=range(mu-4*sigmap,mu+4*sigmap),ylim=range(0,(1+0.1)/(2*sigmap)),type="l",main=title,col=1,xlab="x",ylab="f(x)")
        if(n>0) for(i in 1:n){lines(spline(xfit,dnormp(xfit,mu=mu,sigmap=sigmap,p=p[i])),col=i+1)}
          leg<-paste("p=",p)
        leg<-c(leg,expression(paste("p-> ",infinity)))
        legend(mu+2.5*sigmap, (1-0.1)/(2*sigmap), leg,col=c(1:n+1,1), lty=1)
}
else {
    plot(xfit,dnormp(xfit,mu=mu,sigmap=sigmap,p=p[1]),xlim=range(mu-4*sigmap,mu+4*sigmap),ylim=range(0,(1+0.1)/(2*sigmap)),type="l",main=title,col=1,xlab="x",ylab="f(x)")
        if(n>1) for(i in 2:n){lines(spline(xfit,dnormp(xfit,mu=mu,sigmap=sigmap,p=p[i])),col=i)}
    leg<-paste("p=",p)
    legend(mu+2.5*sigmap, (1-0.1)/(2*sigmap), leg,col=1:n, lty=1)
    }
} 
}


simul.mp<-function(n,m,mu=0,sigmap=1,p=2){
ris<-list()
a<-rnormp(n*m,mu=mu,sigmap=sigmap,p=p)
 a<-matrix(a,n,m)
 xc<-split(a,col(a))
  mat<-matrix(nrow=m,ncol=6)
   for (i in 1:m){
    A<-paramp(xc[[i]])
    for (j in 1:6) mat[i,j]<-A[[j]]
}
ma<-matrix(nrow=2,ncol=5,dimnames=list(c("Mean","Variance"),c("Mean","Mp","Sd","Sp","p")))
  for(j in 1:5){ma[1,j]<-mean(mat[,j]);ma[2,j]<-(length(mat[,j])-1)*var(mat[,j])/length(mat[,j])}
ris$dat<-mat
ris$table<-ma
ris$iter<-sum(mat[,6])
class(ris)<-"simul.mp"
ris
}


print.simul.mp<-function(x,...){
print(x$table)
cat("\nNumber of samples with a difficult convergence:",x$iter,"\n")
invisible(x)
}

plot.simul.mp<-function(x,...){
       mat<-x$dat
       hist(mat[,1],breaks=10,freq=FALSE,xlab="",col="yellow",main=paste("Histogram of",colnames(mat)[1]))  
       op<-par(ask=TRUE)
       for (i in 2:4) {hist(mat[,i],breaks=10,freq=FALSE,xlab ="",col=i,main=paste("Histogram of",colnames(mat)[i]))}
    pval<-mat[,5]
    pval<-pval[pval<=10]
     hist(pval,breaks=seq(1,11,by=1),freq=FALSE,xlab ="",col=5,main=("Histogram of p"))
    par(op)
      }


simul.lmp<-function(n, m, q, data, int=0, sigmap=1, p=2, lp=FALSE){
if (!is.numeric(n)||!is.numeric(m)||!is.numeric(q)||!is.numeric(data)||!is.numeric(int)
||!is.numeric(sigmap)||!is.numeric(p)) 
stop("Non-numeric argument to mathematical function")
if (length(data)!=q) stop("Length of data vector must be q")
name<-paste("x",1:q,sep="")
ris<-matrix(nrow=m,ncol=q+3,dimnames=list(c(1:m),c("(intercept)",name,"Sp","p")))
conv<-vector(length=m)
fr<-"y~"
for(k in 1:q) {fr<-paste(fr,name[k],sep="+")}
frm<-as.formula(fr)
if (lp==FALSE) pp<-NULL
if (lp==TRUE)  pp<-p
for (i in 1:m){
e<-rnormp(n,mu=0,sigmap=sigmap,p=p)
X<-matrix(runif(q*n),nrow=n,ncol=q,dimnames=(list(c(1:n),name)))
y<- int + data %*% t(X) + e
y<-as.vector(y)
X<-as.data.frame(X)
ll<-lmp(formula=frm,data=as.list(X),p=pp)
ll<-summary(ll)
ris[i,]<-c(ll$coeff,ll$sigma,ll$p)
if(lp==FALSE) conv[i]<-c(ll$iter)
}
mat<-ris
  ma<-matrix(nrow=2,ncol=ncol(ris),dimnames=list(c("Mean","Variance"),colnames(ris)))
  for(j in 1:ncol(ma))  {
  ma[1,j]<-mean(mat[,j])
  ma[2,j]<-(m-1)*var(mat[,j])/m
  }
ret<-list(dat=ris,table=ma)
ret$p<-p
ret$par<-c(int,data)
names(ret$par)<-colnames(ma)[1:length(ret$par)]
ret$frm<-frm
if(lp==FALSE) ret$iter<-sum(conv)
ret$lp<-lp
class(ret)<-c("simul.lmp","simul.mp")
ret
}

summary.simul.lmp<-function(object,...){
ans<-object
szca<-nrow(ans$dat)
ans$szca<-szca
class(ans)<-"summary.simul.lmp"
ans
}
  
print.summary.simul.lmp<-function(x,...){
a<-x
cat("Results:\n")
print(a$table)
cat("\nCoefficients:\n")
print(a$par)
cat("\nFormula:\n")
print(a$frm)
cat("\nNumber of samples:\n")
print(as.symbol(a$szca))
cat("\nValue of p:\n")
print(as.symbol(a$p))
if(x$lp==FALSE){
cat("\nNumber of samples with problems on convergence\n")
print(as.symbol(a$iter))
}
cat("\n")
invisible(a)
}

plot.simul.lmp<-function(x,...){
         mat<-x$dat
         hist(mat[,1],breaks=10,freq=FALSE,xlab="",col="blue",main=paste("Histogram of",colnames(mat)[1]))  
         op<-par(ask=TRUE)
         for (i in 2:(ncol(mat)-1)) hist(mat[,i],breaks=10,freq=FALSE,xlab ="",col=rgb(runif(1),runif(1),runif(1)),main=paste("Histogram of",colnames(mat)[i]))
         if(x$lp==FALSE){
              pval<-mat[,ncol(mat)]
              pval<-pval[pval<=10]
              hist(pval,breaks=seq(1,11,by=1),freq=FALSE,xlab ="",col=ncol(mat),main=("Histogram of p"))
         }
         par(op)
      }
