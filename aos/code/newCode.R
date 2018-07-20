myDataGen <- function(n,omega,xi,myV){
    myOut <- rep(0,n)
    X <- rnorm(n,0,1)
    Y <- rnorm(n,xi,myV)
    tmp <- rbinom(n,1, omega)
    return(Y*tmp+X*(1-tmp))
}


myLogLikelihood <- function(myData){
    function(omega,xi,myV){
        mySum <- 0
        for(i in 1:length(myData)){
            mySum = mySum+log((1-omega)*dnorm(myData[i])+omega*dnorm(myData[i],xi,myV))
        }
        return(mySum)
    }
}

myD <- myDataGen(50,1/2,1,0.1)
myF <- myLogLikelihood(myD)
mySmallF <- function(myV){
    myF(1/2, myD[1], myV)
}

tmpX <- seq(0,100,by=0.5)
tmpY <- mySmallF(exp(-tmpX/2))
#plot(exp(-tmpX/2),tmpY)
#plot(tmpX,tmpY)



library(ggplot2)
library(latex2exp)
myP <- ggplot(data.frame(x=tmpX,y=tmpY), aes(x = x, y = y)) +
    geom_line() +
    geom_vline(xintercept=-2*log(0.1),linetype=2)+
    ylab(TeX("Likelihood"))+
    xlab(TeX("$-\\log (\\sigma^2)$"))+
    #annotate("text",x=-2*log(0.1),y=0,label="1[2]",parse=TRUE)+
    annotation_custom("a[b]",xmin=-2*log(0.1),,xmax=-2*log(0.1),,ymin=-60,ymax=-60)+
    theme_classic()
myP
ggsave("myP.png", myP)

