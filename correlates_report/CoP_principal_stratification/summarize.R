source(file="~/Pepe/RA/Functions/FunctionCall.R")

ps.seq<-c('EIA.log10d14overd0','PCA.log10d14overd0')
y.seq<-c('y1','y2')


index1=2
index2=2

ps<-ps.seq[index1]
y.v<-y.seq[index2]




out.VE1.seq<-out.VE2.seq<-out.VE3.seq<-out.VE4.seq<-NULL
out.beta1.seq<-out.beta2.seq<-out.beta3.seq<-out.beta4.seq<-NULL


for (index3 in 1:10){
load(file=paste0("~/RSVcorrelatesAnalysis/Result/outbootCohort1Full_",ps,"_",y.v,"_",index3,".Rdata"))
#print(length(out.VE3))
print(summary(out.ll))
#out.VE1.seq<-rbind(out.VE1.seq,matrix(out.VE1,byrow=T,ncol=length(Su)))
#out.VE2.seq<-rbind(out.VE2.seq,matrix(out.VE2,byrow=T,ncol=length(Su)))
#out.VE3.seq<-try.warning(rbind(out.VE3.seq,matrix(out.VE3,byrow=T,ncol=length(Su))))
#if (inherits(out.VE3.seq,'try-error')) print(index3)
#out.VE4.seq<-rbind(out.VE4.seq,matrix(out.VE4,byrow=T,ncol=length(Su)))
out.VE1.seq<-rbind(out.VE1.seq,out.VE1)
out.VE2.seq<-rbind(out.VE2.seq,out.VE2)
out.VE3.seq<-rbind(out.VE3.seq,out.VE3)
out.VE4.seq<-rbind(out.VE4.seq,out.VE4)
out.beta1.seq<-rbind(out.beta1.seq,out.beta1)
out.beta2.seq<-rbind(out.beta2.seq,out.beta2)
out.beta3.seq<-rbind(out.beta3.seq,out.beta3)
out.beta4.seq<-rbind(out.beta4.seq,out.beta4)
}

save(Su,out.beta1.seq,out.beta2.seq,out.beta3.seq,out.beta4.seq,
out.VE1.seq,out.VE2.seq,out.VE3.seq,out.VE4.seq,file=paste0("~/RSVcorrelatesAnalysis/Result/outbootCohort1Full_",ps,"_",y.v,".Rdata"))



#############
ps.seq<-c('EIA.log10d14overd0','PCA.log10d14overd0')
y.seq<-c('y1','y2')


index1=2
index2=2

ps<-ps.seq[index1]
y.v<-y.seq[index2]




out.VE1.seq<-out.VE2.seq<-out.VE3.seq<-out.VE4.seq<-NULL
out.beta1.seq<-out.beta2.seq<-out.beta3.seq<-out.beta4.seq<-NULL

for (index3 in 1:10){
load(file=paste0("~/RSVcorrelatesAnalysis/Result/outbootCohort1_",ps,"_",y.v,"_",index3,".Rdata"))

#out.VE1.seq<-rbind(out.VE1.seq,matrix(out.VE1,byrow=T,ncol=length(Su)))
#out.VE2.seq<-rbind(out.VE2.seq,matrix(out.VE2,byrow=T,ncol=length(Su)))
#out.VE3.seq<-rbind(out.VE3.seq,matrix(out.VE3,byrow=T,ncol=length(Su)))
#out.VE4.seq<-rbind(out.VE4.seq,matrix(out.VE4,byrow=T,ncol=length(Su)))
out.VE1.seq<-rbind(out.VE1.seq,out.VE1)
out.VE2.seq<-rbind(out.VE2.seq,out.VE2)
out.VE3.seq<-rbind(out.VE3.seq,out.VE3)
out.VE4.seq<-rbind(out.VE4.seq,out.VE4)
out.beta1.seq<-rbind(out.beta1.seq,out.beta1)
out.beta2.seq<-rbind(out.beta2.seq,out.beta2)
out.beta3.seq<-rbind(out.beta3.seq,out.beta3)
out.beta4.seq<-rbind(out.beta4.seq,out.beta4)
}

save(Su,out.beta1.seq,out.beta2.seq,out.beta3.seq,out.beta4.seq,out.VE1.seq,out.VE2.seq,out.VE3.seq,out.VE4.seq,file=paste0("~/RSVcorrelatesAnalysis/Result/outbootCohort1_",ps,"_",y.v,".Rdata"))

#############



ps.seq<-c('EIA.log10d14overd0','PCA.log10d14overd0','RSVA.log10d14','RSVB.log10d14')
bsm.seq<-c('EIA.log10d0','PCA.log10d0','RSVA.log10d0','RSVB.log10d0')
y.seq<-c('y1','y2')


index1=4
index2=2
ps<-ps.seq[index1]
y.v<-y.seq[index2]


out.VE1.seq<-out.VE2.seq<-out.VE3.seq<-out.VE4.seq<-NULL
out.beta1.seq<-out.beta2.seq<-out.beta3.seq<-out.beta4.seq<-NULL

for (index3 in 1:10){
load(file=paste0("~/RSVcorrelatesAnalysis/Result/outbootCC_",ps,"_",y.v,"_",index3,".Rdata"))

#out.VE1.seq<-rbind(out.VE1.seq,matrix(out.VE1,byrow=T,ncol=length(Su)))
#out.VE2.seq<-rbind(out.VE2.seq,matrix(out.VE2,byrow=T,ncol=length(Su)))
#out.VE3.seq<-rbind(out.VE3.seq,matrix(out.VE3,byrow=T,ncol=length(Su)))
#out.VE4.seq<-rbind(out.VE4.seq,matrix(out.VE4,byrow=T,ncol=length(Su)))
out.VE1.seq<-rbind(out.VE1.seq,out.VE1)
out.VE2.seq<-rbind(out.VE2.seq,out.VE2)
out.VE3.seq<-rbind(out.VE3.seq,out.VE3)
out.VE4.seq<-rbind(out.VE4.seq,out.VE4)
out.beta1.seq<-rbind(out.beta1.seq,out.beta1)
out.beta2.seq<-rbind(out.beta2.seq,out.beta2)
out.beta3.seq<-rbind(out.beta3.seq,out.beta3)
out.beta4.seq<-rbind(out.beta4.seq,out.beta4)
}

save(Su,out.beta1.seq,out.beta2.seq,out.beta3.seq,out.beta4.seq,
out.VE1.seq,out.VE2.seq,out.VE3.seq,out.VE4.seq,out.VE1.seq,out.VE2.seq,out.VE3.seq,out.VE4.seq,file=paste0("~/RSVcorrelatesAnalysis/Result/outbootCC_",ps,"_",y.v,".Rdata"))
