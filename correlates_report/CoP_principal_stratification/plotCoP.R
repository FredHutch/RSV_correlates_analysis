setwd("C:/All_Files/1yingsstuff/RSVcorrelatesAnalysis/Result")
library(kyotil)

index1=2;index2=2

mypostscript(file=paste0("../Figure/fig",ps.seq[index1],'_',y.seq[index2]),width=6,height=6)


ps.seq<-c('EIA.log10d14overd0','PCA.log10d14overd0')
y.seq<-c('y1','y2')

ps.name<-c("EIA d14-d0","PCA d14-d0")
y.name<-c("endpoint 1", "endpoint 2")
ps<-ps.seq[index1]
y.v<-y.seq[index2]

load(file=paste0("outCohort_",ps,"_",y.v,".Rdata"))


plot(VE1$Su,VE1$VE,type='l',ylim=c(0.1,0.9),xlab=quote(s[1]),ylab=quote(VE(s[1])),main=paste0(ps.name[index1],',',y.name[index2]))
lines(VE2$Su,VE2$VE,col=2)
lines(VE3$Su,VE3$VE,col=3)
lines(VE4$Su,VE4$VE,col=4)


ps.seq<-c('EIA.log10d14overd0','PCA.log10d14overd0','RSVA.log10d14','RSVB.log10d14')
bsm.seq<-c('EIA.log10d0','PCA.log10d0','RSVA.log10d0','RSVB.log10d0')
y.seq<-c('y1','y2')

ps<-ps.seq[index1]
y.v<-y.seq[index2]


load(file=paste0("outCC_",ps,"_",y.v,".Rdata"))

lines(VE1$Su,VE1$VE,col=1,lty=2)
#lines(VE2$Su,VE2$VE,col=2,lty=2)
lines(VE3$Su,VE3$VE,col=3,lty=2)
#lines(VE4$Su,VE4$VE,col=4,lty=2)

legend("bottomright",col=c(1,3,1,3),lty=c(1,1,2,2),c("Cohort, w/o cov-adj","Cohort, cov-adj","Case-Control, w/o cov-adj","Case-Control, cov-adj"),bty='n',cex=0.8)
dev.off()


index1=3;index2=2

mypostscript(file=paste0("../Figure/fig",ps.seq[index1],'_',y.seq[index2]),width=6,height=6)

ps.seq<-c('EIA.log10d14overd0','PCA.log10d14overd0','RSVA.log10d14','RSVB.log10d14')
bsm.seq<-c('EIA.log10d0','PCA.log10d0','RSVA.log10d0','RSVB.log10d0')
y.seq<-c('y1','y2')

ps.name<-c("EIA d14-d0","PCA d14-d0","RSVA d14","RSVB d14")
y.name<-c("endpoint 1", "endpoint 2")

ps<-ps.seq[index1]
y.v<-y.seq[index2]


load(file=paste0("outCC_",ps,"_",y.v,".Rdata"))

plot(VE1$Su,VE1$VE,col=1,type='l',ylim=c(-0.2,1),xlab=quote(s[1]),ylab=quote(VE(s[1])),main=paste0(ps.name[index1],',',y.name[index2]))
#lines(VE2$Su,VE2$VE,col=2,lty=2)
lines(VE3$Su,VE3$VE,col=3,lty=2)
#lines(VE4$Su,VE4$VE,col=4,lty=2)
legend("bottomright",col=c(1,3),lty=c(2,2),c("Case-Control, w/o cov-adj","Case-Control, cov-adj"),bty='n',cex=0.8)
dev.off()
