library(ggplot2)
library(grid)
myResult=read.csv("myNormalResult.csv")
jjj=list()
for(k in 1:9){
    hh=myResult[,k]
    hh2=sort(hh)
    xxx=numeric()
    for(i in 1:length(hh)){
        xxx[i]=qchisq((i-0.5)/length(hh),df=2)
    }
    myD=data.frame(x=xxx,y=hh2)
    jjj[[k]]=ggplot(data=myD)+
        geom_point(aes(x=x,y=y),colour="darkblue")+
        geom_abline(slope=1,intercept=0)+
        theme_bw()
    #plot(xxx,hh2)
    #abline(a=0,b=1)
}

pdf("figure/test.pdf",width=8,height=6)
grid.newpage()
pushViewport(viewport(layout=grid.layout(3,3)))
vplayout=function(x,y){
    viewport(layout.pos.row=x,layout.pos.col=y)
}
for(i in 1:9){
    print(jjj[[i]],vp=vplayout(ceiling(i/3),i-(ceiling(i/3)-1)*3))
}
dev.off()


temp1=strsplit(names(myResult),split="alpha.")
alpha=numeric()
for(i in 1:9){
    alpha[i]=temp1[[i]][2]
}
alpha=as.numeric(alpha)
myN=NULL
for(i in 1:9){
    myN[i]=strsplit(temp1[[i]][1],split="n.")[[1]][[2]]  
}
myN=as.numeric(myN)

ama=data.frame(x=NULL,alpha=NULL,n=NULL)
for(i in 1:9){
    temp1=sort(myResult[,i])
    xxx=NULL
    for(j in 1:length(temp1)){
        xxx[j]=qchisq((j-0.5)/length(temp1),df=2)
    }
    temp=data.frame(x=xxx,y=temp1,alpha=alpha[i],n=myN[i])
    ama=rbind(ama,temp)
}

myQQPlot=ggplot(data=ama)+geom_point(mapping=aes(x=x,y=y),colour="darkblue")+
    geom_abline(slope = 1,intercept = 0)+
    facet_grid(n~alpha,labeller=label_bquote(cols=alpha== .(alpha),rows=n== .(n)))+
    labs(title="Q-Q plot",x="Theoretical Quantile",y="Sample Quantile")+
    theme_minimal()+
    theme(panel.background=element_rect(colour="gray20"),
          axis.ticks.length=unit(0.1,"cm"),
          axis.ticks=element_line(colour="black",size=2,linetype="dotted",lineend="round"))
#?plotmath
pdf("../figure/myQQPlotNormal.pdf",width=4.5,height=4)
#grid.newpage()
print(myQQPlot)
dev.off()

hist(myResult[,9])
myPvalue=pchisq(myResult[,9],df=2)
hist(myPvalue)
sum(myPvalue>0.95)/length(myResult[,1])
