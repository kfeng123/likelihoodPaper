n=60
p=120

Bai_test_pvalue=function(X){
    S=cov(X)
    X_bar=colMeans(X)
    dim(X_bar)=c(p,1)
    T=(n*t(X_bar)%*%X_bar-sum(diag(S)))/(2*(n-1)*n/(n-2)/(n+1)*(sum(diag(S%*%S))-sum(diag(S))^2/(n-1)))^(1/2)
    1-pnorm(T)
}
new_likelihood_pvalue=function(X){
    B=100
    calc_T=function(X){
        S=cov(X)
        X_bar=colMeans(X)
        dim(X_bar)=c(p,1)
        t(X_bar)%*%solve(diag(p)/n+S)%*%X_bar
    }
    myT=calc_T(X)
    random_stat=NULL
    for(i in 1:B){
        the_sign=sample(c(-1,1),n,replace=TRUE) 
        the_sign=diag(the_sign)
        random_stat[i]=calc_T(the_sign%*%X)
    }
    sum(random_stat>as.numeric(myT))/B
}

jjie1=NULL
jjie2=NULL
k=0
for(c in seq(0,1,0.05)){
    k=k+1
    hh=NULL
    newHH=NULL
    for(i in 1:1000){
        X=rnorm(n*p,0,1)
        dim(X)=c(n,p)
        X[,1]=X[,1]/sqrt(100)
        X[,1]=X[,1]+1*c
        X[,2:p]=X[,2:p]+0.01*c
        hh[i]=Bai_test_pvalue(X)
        newHH[i]=new_likelihood_pvalue(X)
    }
    #hist(hh)
    #hist(newHH)
    jjie1[k]=sum(hh<0.05)/1000
    jjie2[k]=sum(newHH<0.05)/1000
}
myD=data.frame(bai=jjie1,new=jjie2)
write.csv(myD,"HDPower2.csv",row.names=FALSE)
