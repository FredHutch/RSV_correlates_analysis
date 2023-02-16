rm(list=ls())    
save.results.to="input"    

index1<-as.integer(commandArgs(trailingOnly=T)[1])
index2<-as.integer(commandArgs(trailingOnly=T)[2])



#index1=3
#index2=1


library(kyotil)
#library(osDesign)
#library(survey)
library(RSVcorr)
library(splines)
library(nnet)
library(dummies)
# options for svyglm to deal with strata with counts of 1
options(survey.lonely.psu = "adjust")
#
dat.wide.v$riskScore.mat.endpoint1=scale(dat.wide.v$riskScore.mat.endpoint1)
dat.wide.v$riskScore.mat.endpoint2=scale(dat.wide.v$riskScore.mat.endpoint2)
dat.wide$riskScore.mat.endpoint1=scale(dat.wide$riskScore.mat.endpoint1)
dat.wide$riskScore.mat.endpoint2=scale(dat.wide$riskScore.mat.endpoint2)


#cohort.seq=c("ppimm","ppimi")
#cohort=cohort.seq[index0]
## repeat the analysis for the two PP-IMM cohorts
##for (cohort in c("PP-IMM-M", "PP-IMM-I")){
#  if (cohort=="ppimm"){
#    dat.wide.ppim <- subset(dat.wide, ppimmfl=="Y")
#  } else {
#    dat.wide.ppim <- subset(dat.wide, ppimifl=="Y")
#  }

### use ITT population
dat.wide.ppim<-dat.wide
 
#> nrow(dat.wide)
#[1] 2337
#> nrow(dat.wide.v)
#[1] 417
#
#> table(dat.wide$ppimmfl)
#
#   N    Y 
# 198 2139 
#> table(dat.wide.v$ppimmfl)
#
#  N   Y 
# 43 374 

#> table(dat.wide$vacc2birthLESSTHAN30)
#
#   0    1 
#2059  278 
#> table(dat.wide.ppim$riskScoreBin.mat.endpoint1)
#
#   0    1 
#1138 1001 

### Need to create a new sampled variable that map phase-2 data

dat.wide.ppim$sampled<-ifelse(dat.wide.ppim$pair.id%in%dat.wide.v$pair.id,'Y','N')

tmp = factor(sub("_case", "", sub("_ctrl",
        "", sub("_vacc","", sub("_plac","", dat.wide.ppim$stratum)))))
    dat.wide.ppim$group = as.integer(tmp)
   

source(file="~/RSVcorrelatesAnalysis/YingFunctions/NonparFun_NoW.R")
source(file="~/RSVcorrelatesAnalysis/YingFunctions/Fun_2020.R")
source(file="~/RSVcorrelatesAnalysis/YingFunctions/function10_15_Full.R")
source(file="~/Pepe/RA/Functions/FunctionCall.R")

 
## endpoints 1 and 2
#  for (y in c("y1", "y2")){
#    # EIA d14 and fold-change, and PCA d14 and fold-change (all on log10 scale)
#    for (ps in c("EIA.log10d14", "EIA.log10d14overd0", "PCA.log10d14", "PCA.log10d14overd0")){
#      
#      formula <- as.formula(paste0(y, " ~ ", ps, " + factor(vacc2birthLESSTHAN30) + factor(riskScoreBin.mat.endpoint", substr(y, 2, 2), ")"))
#      print(formula)
#      
#      riskCurve(formula, bsm=paste0(substr(ps, 1, 9), "d0"), tx="trt", data=dat.wide.ppim, weights="wt", 
#                saveFile=paste0("VEcurve_", cohort, "_", y, "_", ps, ".RData"), saveDir=outDir)
#    }
#  }
#}

ps.seq<-c('EIA.log10d14overd0','PCA.log10d14overd0','RSVA.log10d14','RSVB.log10d14')
bsm.seq<-c('EIA.log10d0','PCA.log10d0','RSVA.log10d0','RSVB.log10d0')
y.seq<-c('y1','y2')

ps<-ps.seq[index1]
y.v<-y.seq[index2]
c.v<-paste0("riskScoreBin.mat.endpoint", substr(y.v, 2, 2))
bsm=bsm.seq[index1]

 # dichotomize risk score to create a small number of baseline strata
  med <- quantile(dat.wide.ppim[,paste0("riskScore.mat.endpoint",substr(y.v,2,2))], probs=0.5)
  dat.wide.ppim[,c.v] <- ifelse(dat.wide.ppim[,paste0("riskScore.mat.endpoint",substr(y.v,2,2))]<=med, 0, 1)
  
############
#### process phase-II data for VE curve computation

 #if (cohort=="ppimm"){
#    dat.wide.v.ppim <- subset(dat.wide.v, ppimmfl=="Y")
#  } else {
#    dat.wide.v.ppim <- subset(dat.wide.v, ppimifl=="Y")
#  }
### use ITT population
dat.wide.v.ppim<-dat.wide.v
dat.wide.v.ppim$S1<-dat.wide.v.ppim[,ps]
dat.wide.v.ppim$Z<-dat.wide.v.ppim$trt
dat.wide.v.ppim$Y<-dat.wide.v.ppim[,y.v]



qq<-quantile(dat.wide.v.ppim[,bsm],c(.25,.5,.75),na.rm=T)

dat.wide.v.ppim$W<-as.numeric(cut(dat.wide.v.ppim[,bsm],c(-Inf,qq,Inf)))
dat.wide.v.ppim$S1[dat.wide.v.ppim$Z==0]<-NA
dat.wide.v.ppim$S1[is.na(dat.wide.v.ppim$W)]<-NA


dat.wide.v.ppim$group = factor(sub("_case", "", sub("_ctrl",
        "", sub("_vacc","", sub("_plac","", dat.wide.v.ppim$stratum)))))
    dat.wide.v.ppim$X = as.integer(dat.wide.v.ppim$group)
 
dat.wide.v.ppim$VB=dat.wide.v.ppim$vacc2birthLESSTHAN30  
dat.wide.v.ppim$RS=ifelse(dat.wide.v.ppim[,paste0("riskScore.mat.endpoint",substr(y.v,2,2))]<=med, 0, 1)
dat.wide.v.ppim$weights=dat.wide.v.ppim$wt  

######  
  
dat.wide.ppim$VB=dat.wide.ppim$vacc2birthLESSTHAN30  
dat.wide.ppim$RS=dat.wide.ppim[,c.v]



dat.wide.ppim<-merge(dat.wide.ppim,dat.wide.v.ppim[,c('pair.id','W','S1')],by="pair.id",all.x=T)

#S1<-dat.wide.ppim[,ps]
S1<-dat.wide.ppim$S1
Z<-dat.wide.ppim$trt
Y<-dat.wide.ppim[,y.v]

#W<-as.numeric(cut(dat.wide.ppim[,bsm],c(-Inf,qq,Inf)))
#Wu<-unique(sort(W))

W<-dat.wide.ppim$W
Wu<-unique(sort(W))

#dat.wide.ppim$W=W
#datin<-data.frame(S1=S1,Z=Z,Y=Y,X=X,W=W)

#fitR=glm(Y~Z+S1+I(Z*S1)+I(W==2)+I(W==3)+I(W==4),family=binomial(link=probit),data=dat.wide.ppim)
#beta=coef(fitR)


#S1[Z==0]<-NA
#S1[is.na(W)]<-NA
#S1[dat.wide.ppim$sampled=='N']<-NA  #### excluding markers not included in case-control sample

Su=sort(unique(S1))
delta1=as.numeric(!is.na(W))
delta2=as.numeric(!is.na(S1))

X<-dat.wide.ppim$group
Xu<-unique(sort(X))

Sout.NC.All.X<-covDFcc.Nonpar2.X.New(S1[delta2==1],W[delta2==1],X[delta2==1],Wu=Wu,Xu=Xu,weights=rep(1,sum(delta2)))

WoutA.NC.All.X<-covDFcc.Nonpar2.Full(W[delta2==1],X[delta2==1],Wu=Xu)

source(file="~/RSVcorrelatesAnalysis/YingFunctions/functionPar2020.R")

 # fit1<-try.error(
#   EM.cc.CPV.Probit.Small.Nonpar2.NC.SubWA.X.Short.New(Z,Sout.NC.All.X,WoutA.NC.All.X,S1,W,Y,X,Wu,Xu,delta1,delta2,pd2.YZXW,c(0,0,0,0)))
#
#
#W<-dat.wide.ppim$vacc2birthLESSTHAN30*100+dat.wide.ppim$riskScoreBin.mat.endpoint1*10+as.numeric(cut(dat.wide.ppim[,bsm],c(-Inf,qq,Inf)))
#Wu<-unique(sort(W))
#
#
#datin<-data.frame(S1=S1,Z=Z,Y=Y,X=X,W=W)
#
#Sout.NC.All.X<-covDFcc.Nonpar2.X.New(S1[delta2==1],W[delta2==1],X[delta2==1],Wu=Wu,Xu=Xu,weights=rep(1,sum(delta2)))
#
#WoutA.NC.All.X<-covDFcc.Nonpar2.Full(W[delta2==1],X[delta2==1],Wu=Xu)
#
#source(file="~/RSVcorrelatesAnalysis/YingFunctions/functionPar.R")
#
#
#
#Sout=Sout.NC.All.X;Wout=WoutA.NC.All.X;varlist=c("VB","RS");beta=rep(0,6)
#Sout=Sout.NC.All.X;Wout=WoutA.NC.All.X;varlist=NULL;beta=rep(0,4)



#datin<-data.frame(S1=S1,Z=Z,Y=Y,X=X,W=W,weights=dat.wide.ppim$wt.ppimmfl,VB=dat.wide.ppim$vacc2birthLESSTHAN30,RS=dat.wide.ppim$riskScoreBin.mat.endpoint1)
#datin.nomis<-datin[!is.na(datin$S1),]




### No covariate adjustment, baseline marker not affecting risk 
fit1<-try.error(
   EM.cc.CPV.Probit.Small.Nonpar2.NC.SubWA.X.Short.New.2020.X(Z,Sout.NC.All.X,WoutA.NC.All.X,S1,W,Y,X,Wu,Xu,delta1,delta2,pd2.YZXW,beta=c(0,0,0,0)))


VE1<-integVE(Su,Wu=Wu,dat.wide.v.ppim,varlist=NULL,beta=fit1)

### No covariate adjustment, allow baseline marker to affect risk
fit2<-try.error(
   EM.cc.CPV.Probit.Small.Nonpar2.NC.SubWA.X.Short.New.2020.X(Z,Sout.NC.All.X,WoutA.NC.All.X,S1,W,Y,X,Wu,Xu,delta1,delta2,pd2.YZXW,varlist=c("W2","W3","W4"),beta=rep(0,7)))

VE2<-integVE(Su,Wu=Wu,dat.wide.v.ppim,varlist=c("W2","W3","W4"),beta=fit2)


### With covariate adjustment, baseline marker not affecting risk
## expand X
X<-dat.wide.ppim$VB*1000+dat.wide.ppim$RS*100+dat.wide.ppim$group
Xu<-unique(sort(X))

Sout.NC.All.X<-covDFcc.Nonpar2.X.New(S1[delta2==1],W[delta2==1],X[delta2==1],Wu=Wu,Xu=Xu,weights=rep(1,sum(delta2)))

WoutA.NC.All.X<-covDFcc.Nonpar2.Full(W[delta2==1],X[delta2==1],Wu=Xu)

source(file="~/RSVcorrelatesAnalysis/YingFunctions/functionPar2020.R")

fit3<-try.error(
   EM.cc.CPV.Probit.Small.Nonpar2.NC.SubWA.X.Short.New.2020.X(Z,Sout.NC.All.X,WoutA.NC.All.X,S1,W,Y,X,Wu,Xu,delta1,delta2,pd2.YZXW,varlist=c("VB","RS"),beta=rep(0,6)))


WX<-dat.wide.ppim$VB*100+dat.wide.ppim$RS*10+dat.wide.ppim$W
WXu<-unique(sort(WX))

dat.wide.v.ppim$W=dat.wide.v.ppim$VB*100+dat.wide.v.ppim$RS*10+dat.wide.v.ppim$W

VE3<-integVE(Su,Wu=WXu,dat.wide.v.ppim,varlist=c("VB","RS"),beta=fit3)

### With covariate adjustment, baseline marker affecting risk

fit4<-try.error(
   EM.cc.CPV.Probit.Small.Nonpar2.NC.SubWA.X.Short.New.2020.X(Z,Sout.NC.All.X,WoutA.NC.All.X,S1,W,Y,X,Wu,Xu,delta1,delta2,pd2.YZXW,varlist=c("W2","W3","W4","VB","RS"),
   beta=rep(0,9)))


VE4<-integVE(Su,Wu=WXu,dat.wide.v.ppim,varlist=c("W2","W3","W4","VB","RS"),beta=fit4)


save(fit1,fit2,fit3,fit4,VE1,VE2,VE3,VE4,file=paste0("~/RSVcorrelatesAnalysis/Result/outCC_",ps,"_",y.v,".Rdata"))
