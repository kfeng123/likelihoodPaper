library(cubature)
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
        myScale=adaptIntegrate(function(x){
            myPrior(x[1],x[2])*myLikelihood(x[1],x[2])
        },lowerLimit=c(-2,0),upperLimit = c(2,5))
        
        myPosterior=function(myMu,mySigma){
            myPrior(myMu,mySigma)*myLikelihood(myMu,mySigma)/myScale$integral
        }
        
        myStatistic=adaptIntegrate(function(x){
            myPosterior(x[1],x[2])*myLikelihood(x[1],x[2])
        },lowerLimit=c(-2,0),upperLimit = c(2,5))$integral
        
        hh[k]=2*log(myStatistic)+2*log(2)
    
        hh[k]=max(hh[k],0)
    }
    
    jj=paste("n=",n,"alpha=",alpha,sep="")
    myResult[,jj]=hh
} 
write.csv(myResult,file="myResult.csv",row.names=FALSE)
