library(ggplot2)
library(grid)
myD=read.csv("HDPower2.csv")
plot(myD[,1])
plot(myD[,2])
myD$c=seq(0,1,0.05)
myPowerFigure=ggplot(myD)+
    geom_line(aes(x=c,y=bai,colour="Bai Test",linetype="Bai Test"))+
    geom_line(mapping=aes(x=c,y=new,colour="New Test",linetype="New Test"))+
    scale_colour_manual(name="",values=c("Bai Test"="black","New Test"="blue"))+
    scale_linetype_manual(name="",values=c("Bai Test"="dashed","New Test"="solid"))+
    labs(title="Power Comparison",x="c",y="Power")+
    theme_bw()
pdf("../figure/myPowerFigure.pdf",width=4.5,height=3.3)
print(myPowerFigure)
dev.off()
