
library(cubature)
library(stats4)
library(mvtnorm)
# Normal approximation

K=1000
myResult=data.frame(row.names = seq(1,K))

for(n in c(50,100,200))for(alpha in c(0.1,0.5,0.9)){
    tempResult=list()
    #normal mixture
    #parameter
    #alpha=0.5
    #n=50
    mu=0
    sigma=2
    #bounded parameter space
    maxSigma=5
    hh=numeric()
    for(k in 1:K){
        #generate data
        X=numeric()
        for(i in 1:n){
            temp1=rbinom(1,1,alpha)
            if(temp1==0){
                X=c(X,rnorm(1,mu,1))  
            }else{
                X=c(X,rnorm(1,0,sigma)) 
            }
        }
        #full
        #prior
        myPrior=function(myMu,mySigma){
            dnorm(myMu,0,1)*dchisq(mySigma,df=3)#/pchisq(10,df=3)
        }
        myLikelihood=function(myMu,mySigma){
            temp=(1-alpha)*dnorm(X,myMu,1)+alpha*dnorm(X,myMu,mySigma)
            temp2=(1-alpha)*dnorm(X,mu,1)+alpha*dnorm(X,mu,sigma)
            #exp(sum(log(temp)))
            prod(temp/temp2)   
        }
        
        myPosterior=function(myMu,mySigma){
            myPrior(myMu,mySigma)*myLikelihood(myMu,mySigma)
        }
        #temp2=seq(0.0001,5,length.out = 1000)
        #jjj=sapply(temp2,function(x){myPosterior(0,x)})
        #jjj2=sapply(temp2,function(x){myLikelihood(0,x)})
        #plot(temp2,jjj,type="l")
        #plot(temp2,jjj2,type="l")
               
        
        mleTemp=mle(minuslogl = function(myMu,logSigma){-log(myPosterior(myMu,exp(logSigma)))}
            ,start = list(myMu=1,logSigma=log(2))) 
        estMu=mleTemp@coef["myMu"]
        estSigma=exp(mleTemp@coef["logSigma"])
        I11=function(x,mu,v){
            ((1-alpha)*((x-mu)^2-1)*dnorm(x,mu,1)+alpha*((x-mu)^2/v^2-v^{-1})*dnorm(x,mu,v))/
                ((1-alpha)*dnorm(x,mu,1)+alpha*dnorm(x,mu,v))-
                ((
                    (1-alpha)*(x-mu)*dnorm(x,mu,1)+alpha*(x-mu)/v*dnorm(x,mu,v)
                )/((1-alpha)*dnorm(x,mu,1)+alpha*dnorm(x,mu,v)))^2
        }
        I12=function(x,mu,v){
            (
                (3*alpha*(mu-x)/2/v^2-alpha*(mu-x)^3/2/v^3)*dnorm(x,mu,v)
            )/((1-alpha)*dnorm(x,mu,1)+alpha*dnorm(x,mu,v))-
                ((
                    alpha*((mu-x)^2/2/v^2-1/2/v)*dnorm(x,mu,v)*(
                        (1-alpha)*(x-mu)*dnorm(x,mu,1)+alpha*(x-mu)/v*dnorm(x,mu,v)
                    )
                )/((1-alpha)*dnorm(x,mu,1)+alpha*dnorm(x,mu,v)))^2
        }
        I22=function(x,mu,v){
            (
                alpha*(3/4/v^2-3*(x-mu)^2/2/v^3+(x-mu)^4/4/v^4)*dnorm(x,mu,v)
            )/((1-alpha)*dnorm(x,mu,1)+alpha*dnorm(x,mu,v))-
                ((
                    alpha*((mu-x)^2/2/v^2-1/2/v)*dnorm(x,mu,v)
                )/((1-alpha)*dnorm(x,mu,1)+alpha*dnorm(x,mu,v)))^2
        }
        myI11=mean(sapply(X,function(x){I11(x,estMu,estSigma)}))
        myI12=mean(sapply(X,function(x){I12(x,estMu,estSigma)}))
        myI22=mean(sapply(X,function(x){I22(x,estMu,estSigma)}))
        myVar=solve(-matrix(c(myI11,myI12,myI12,myI22),nrow=2))/n
        #myWeight=function(mu){
        #    dmvnorm(mu,mean=c(estMu,estSigma),sigma=myVar)
        #}
        #myStatistic=adaptIntegrate(function(x){
        #    myWeight(x)*myLikelihood(x[1],x[2])
        #},lowerLimit=c(-2,0),upperLimit = c(2,5))$integral
        ttt1=eigen(myVar)$values
        ttt2=eigen(myVar)$vectors
        ttt1[ttt1<0]=1e-6
        myVar=ttt2%*%diag(ttt1)%*%t(ttt2)
        tttmp=rmvnorm(100,mean=c(estMu,estSigma),sigma=myVar)
        tttmp[,2][tttmp[,2]<0]=1e-6
        myStatistic=mean(apply(tttmp,1,function(x){myLikelihood(x[1],x[2])}))
        
         
        hh[k]=2*log(myStatistic)+2*log(2)
    
        hh[k]=max(hh[k],0)
    }
    
    jj=paste("n=",n,"alpha=",alpha,sep="")
    myResult[,jj]=hh
} 
write.csv(myResult,file="myNormalResult.csv",row.names=FALSE)
