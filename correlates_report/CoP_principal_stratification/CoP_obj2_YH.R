rm(list=ls())    
save.results.to="input"    
library(kyotil)
#library(osDesign)
#library(survey)
library(RSVcorr)
library(splines)
library(nnet)
# options for svyglm to deal with strata with counts of 1
options(survey.lonely.psu = "adjust")
#
dat.wide.v$riskScore.mat.endpoint1=scale(dat.wide.v$riskScore.mat.endpoint1)
dat.wide.v$riskScore.mat.endpoint2=scale(dat.wide.v$riskScore.mat.endpoint2)
dat.wide$riskScore.mat.endpoint1=scale(dat.wide$riskScore.mat.endpoint1)
dat.wide$riskScore.mat.endpoint2=scale(dat.wide$riskScore.mat.endpoint2)


cohort="PP-IMM-M"
# repeat the analysis for the two PP-IMM cohorts
#for (cohort in c("PP-IMM-M", "PP-IMM-I")){
  if (cohort=="PP-IMM-M"){
    dat.wide.ppim <- subset(dat.wide, ppimmfl=="Y")
  } else {
    dat.wide.ppim <- subset(dat.wide, ppimifl=="Y")
  }

  # dichotomize risk score to create a small number of baseline strata
  med <- quantile(dat.wide.ppim$riskScore.mat.endpoint1, probs=0.5)
  dat.wide.ppim$riskScoreBin.mat.endpoint1 <- ifelse(dat.wide.ppim$riskScore.mat.endpoint1<=med, 0, 1)
  med <- quantile(dat.wide.ppim$riskScore.mat.endpoint2, probs=0.5)
  dat.wide.ppim$riskScoreBin.mat.endpoint2 <- ifelse(dat.wide.ppim$riskScore.mat.endpoint2<=med, 0, 1)
  
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

ps<-"EIA.log10d14overd0"
bsm=paste0(substr(ps, 1, 9),"d0")
qq<-quantile(dat.wide.ppim[,bsm],c(.25,.5,.75))


S1<-dat.wide.ppim[,ps]
Z<-dat.wide.ppim$trt
Y<-dat.wide.ppim$y1
X<-rep(1,length(Y))
Xu<-unique(X)

W<-as.numeric(cut(dat.wide.ppim[,bsm],c(-Inf,qq,Inf)))
Wu<-unique(sort(W))

datin<-data.frame(S1=S1,Z=Z,Y=Y,X=X,W=W)

#fitR=glm(Y~Z+S1+I(Z*S1)+I(W==2)+I(W==3)+I(W==4),family=binomial(link=probit),data=dat.wide.ppim)
#beta=coef(fitR)


S1[Z==0]<-NA

Su=sort(unique(S1))
delta1=as.numeric(!is.na(W))
delta2=as.numeric(!is.na(S1))


Sout.NC.All.X<-covDFcc.Nonpar2.X.New(S1[delta2==1],W[delta2==1],X[delta2==1],Wu=Wu,Xu=Xu,weights=rep(1,sum(delta2)))

WoutA.NC.All.X<-covDFcc.Nonpar2.Full(W[delta2==1],X[delta2==1],Wu=Xu)

source(file="~/RSVcorrelatesAnalysis/YingFunctions/functionPar.R")

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



datin<-data.frame(S1=S1,Z=Z,Y=Y,X=X,W=W)
datin.nomis<-datin[!is.na(datin$S1),]
datin.nomis$weights=1

### No covariate adjustment, baseline marker not affecting risk 
fit1<-try.error(
   EM.cc.CPV.Probit.Small.Nonpar2.NC.SubWA.X.Short.New.2020(Z,Sout.NC.All.X,WoutA.NC.All.X,S1,W,Y,X,Wu,Xu,delta1,delta2,pd2.YZXW,beta=c(0,0,0,0)))


#fit11<-try.error(
#   EM.cc.CPV.Probit.Small.Nonpar2.NC.SubWA.X.Short.New.2020.X(Z,Sout.NC.All.X,WoutA.NC.All.X,S1,W,Y,X,Wu,Xu,delta1,delta2,pd2.YZXW,beta=c(0,0,0,0)))


VE1<-integVE(Su,Wu,datin.nomis,varlist=NULL,beta=fit1)

### No covariate adjustment, allow baseline marker to affect risk
fit2<-try.error(
   EM.cc.CPV.Probit.Small.Nonpar2.NC.SubWA.X.Short.New.2020(Z,Sout.NC.All.X,WoutA.NC.All.X,S1,W,Y,X,Wu,Xu,delta1,delta2,pd2.YZXW,varlist=c("W2","W3","W4"),beta=rep(0,7)))


datin<-data.frame(S1=S1,Z=Z,Y=Y,X=X,W=W)

VE2<-integVE(Su,Wu,datin.nomis,varlist=c("W2","W3","W4"),beta=fit2)


#### With covariate adjustment, baseline marker not affecting risk 
W<-dat.wide.ppim$vacc2birthLESSTHAN30*100+dat.wide.ppim$riskScoreBin.mat.endpoint1*10+as.numeric(cut(dat.wide.ppim[,bsm],c(-Inf,qq,Inf)))
Wu<-unique(sort(W))
#
Sout.NC.All.X<-covDFcc.Nonpar2.X.New(S1[delta2==1],W[delta2==1],X[delta2==1],Wu=Wu,Xu=Xu,weights=rep(1,sum(delta2)))
#
WoutA.NC.All.X<-covDFcc.Nonpar2.Full(W[delta2==1],X[delta2==1],Wu=Xu)
#
source(file="~/RSVcorrelatesAnalysis/YingFunctions/functionPar.R")
#
#

datin<-data.frame(S1=S1,Z=Z,Y=Y,X=X,W=W)
datin.nomis<-datin[!is.na(datin$S1),]
datin.nomis$weights=1


fit3<-try.error(
   EM.cc.CPV.Probit.Small.Nonpar2.NC.SubWA.X.Short.New.2020(Z,Sout.NC.All.X,WoutA.NC.All.X,S1,W,Y,X,Wu,Xu,delta1,delta2,pd2.YZXW,varlist=c("VB","RS"),beta=rep(0,6)))

VE3<-integVE(Su,Wu,datin.nomis,varlist=c("VB","RS"),beta=fit3)



fit4<-try.error(
   EM.cc.CPV.Probit.Small.Nonpar2.NC.SubWA.X.Short.New.2020(Z,Sout.NC.All.X,WoutA.NC.All.X,S1,W,Y,X,Wu,Xu,delta1,delta2,pd2.YZXW,varlist=c("VB","RS","W2","W3","W4"),beta=rep(0,9)))

VE4<-integVE(Su,Wu,datin.nomis,varlist=c("VB","RS","W2","W3","W4"),beta=fit4)
