# @ by Youyi Fong
# start R in the correlates_report folder or make sure working directory is here
rm(list=ls())    
save.results.to="input"    
library(RSVcorr); stopifnot(packageVersion("RSVcorr")>="2020.10.31")
library(chngpt)
library(osDesign)
library(survey)
# options for svyglm to deal with strata with counts of 1
options(survey.lonely.psu = "adjust")
library(mgcv) # this needs to come before kyotil because it also defines %.%
library(kyotil); stopifnot(packageVersion("kyotil")>="2020.10.12")



#########################################################################################
# baseline CoVE 

var.ind=6
for (y in c("y1","y2","y3")) {
#for (ind in c("a","b")) { # ind a: EIA/PCA use phase 1 data; ind b: EIA/PCA use phase 2 data
#y="y3"; x="RSVA.d0"
    
    markers=c(t(outer(assays,".d0", paste0))); markers; names(markers)=markers    
    fits=list()
    for (x in markers) {
        formula=paste0(y, "~trt*", ifelse(contain(x,"log10"),x,"log10("%.%x%.%")"))
        # add maternal risk scores
        if (y=="y1") formula=paste0(formula,"+riskScore.mat.endpoint1") else formula=paste0(formula,"+riskScore.mat.endpoint2")
        # adj vacc2birth for maternal timepoints
        if(!contain(x,"cord")) formula=paste0(formula,"+vacc2birth")
        formula=as.formula(formula); print(formula)
        #if(ind=="b" | startsWith(x,"RSV")) {
        if(startsWith(x,"RSV")) {
            if (y=="y1") {
                fit=tps.rsv(formula=formula, trt.grp=2)
            } else {
                dstrat<-svydesign(id=~1,strata=~stratum, weights=~wt, data=dat.wide.v)
                fit=svyglm(formula, design=dstrat, family="binomial")
            }                                
            #fit=glm(formula, subset(dat.wide.v, trt==.trt), family="binomial", weights=subset(dat.wide.v, trt==.trt)$wt)
        } else {
            fit=glm(formula=formula, dat.wide, family="binomial")
        }    
        print(summary(fit))
        fits[[x]]=fit
    }
    
    res=c(getFormattedSummary(fits, robust=F, exp=T, rows=var.ind))    
    tab=matrix(res, ncol=length(assays), dimnames=list(paste0(ifelse(contain(times[1],"log10"),"","log10 "),times[1]), assays)); print(tab)
    mytex(tab, file=paste0("tables/CoVE_baseline_", y), align="c", input.folder=save.results.to)
    #mytex(tab, file=paste0("tables/CoVE_baseline", ind, "_", y), align="c", input.folder=save.results.to)
        
}
#}



# including all four baseline markers
var.ind=9
for (y in c("y1","y2","y3")) {
#y="y1"; x="RSVA.d0"
    
    markers=c(t(outer(assays,".d0", paste0))); markers; names(markers)=markers    
    fits=list()
    for (x in markers) {
        formula=paste0(y, "~log10(RSVA.d0)+log10(RSVB.d0)+log10(EIA.d0)+log10(PCA.d0)+trt*", ifelse(contain(x,"log10"),x,"log10("%.%x%.%")"))
        # add maternal risk scores
        if (y=="y1") formula=paste0(formula,"+riskScore.mat.endpoint1") else formula=paste0(formula,"+riskScore.mat.endpoint2")
        # adj vacc2birth for maternal timepoints
        if(!contain(x,"cord")) formula=paste0(formula,"+vacc2birth")
        formula=as.formula(formula); print(formula)
        if (y=="y1") {
            fit=tps.rsv(formula=formula, trt.grp=2)
        } else {
            dstrat<-svydesign(id=~1,strata=~stratum, weights=~wt, data=dat.wide.v)
            fit=svyglm(formula, design=dstrat, family="binomial")
        }                                
        print(summary(fit))
        fits[[x]]=fit
    }
    
    res=c(getFormattedSummary(fits, robust=F, exp=T, rows=var.ind))    
    tab=matrix(res, ncol=length(assays), dimnames=list(paste0(ifelse(contain(times[1],"log10"),"","log10 "),times[1]), assays)); print(tab)
    mytex(tab, file=paste0("tables/CoVE_baseline_4_", y), align="c", input.folder=save.results.to)
    #mytex(tab, file=paste0("tables/CoVE_baseline", ind, "_", y), align="c", input.folder=save.results.to)
        
}


#########################################################################################
# both arms together, interaction with delta or D14

y="y2"
for (ind in 1:2) { 
for (t in c("log10d14overd0","log10d14")) {    
    fits=list()
    for (x in assays) {
        formula=paste0(y, ifelse(ind==1, "~trt*", "~trt+"), x,".",t)    
        if (y=="y1") formula=paste0(formula,"+riskScore.mat.endpoint1") else formula=paste0(formula,"+riskScore.mat.endpoint2") # add maternal risk scores    
        if(!contain(x,"cord")) formula=paste0(formula,"+vacc2birth") # adj vacc2birth for maternal timepoints
        formula=as.formula(formula); print(formula)
        
        if(startsWith(x,"RSV")) {
            if (y=="y1") {
                fit=tps.rsv(formula=formula, trt.grp=2)
            } else {
                dstrat<-svydesign(id=~1,strata=~stratum, weights=~wt, data=dat.wide.v)
                fit=svyglm(formula, design=dstrat, family="binomial")
            }                                
            #fit=glm(formula, subset(dat.wide.v, trt==.trt), family="binomial", weights=subset(dat.wide.v, trt==.trt)$wt)
        } else {
            fit=glm(formula=formula, dat.wide, family="binomial")
        }    
        #print(summary(fit))
        fits[[x]]=fit
    }
    
    tab=getFormattedSummary(fits, robust=F, exp=T)    
    tab=tab[-1,]
    rownames(tab)=sub("RSVA.","",rownames(tab))
    mytex(tab, file=paste0("tables/",ifelse(ind==1, "interaction", "main"),"_", y, "_", t), align="c", input.folder=save.results.to)
}
}



#########################################################################################
# Mediation analysis


# the two formulae are for with and without adjusting baseline markers
formula.effect.1=y2~trt+riskScore.mat.endpoint2+vacc2birth
formula.effect.2=y2~trt+riskScore.mat.endpoint2+vacc2birth+EIA.log10d0
formula.mediators=~EIA.log10d14

myboxplot(EIA.log10d14~trt, dat.wide, log="")
fit.1=iorw(formula.effect.1, formula.mediators, dat.wide, family=binomial(), nboot=0, numCores=1, save.steps=T)
fit.1

fit.2=iorw(formula.effect.2, formula.mediators, dat.wide, family=binomial(), nboot=0, numCores=1, save.steps=T)
fit.2

tab=rbind(fit.1$est, fit.2$est)
rownames(tab)=c("without adjusting baseline marker", "adjusting baseline marker")
tab
mytex(tab, file.name=paste0("mediation"), align="c", input.folder=save.results.to, save2input.only=TRUE)

# make a plot to compare RSV dataset and Cowling et al flu dataset
myfigure (mfrow=c(1,2))
    myboxplot(EIA.log10d14~trt, dat.wide, log="", ylab="EIA log10 D14",   main="RSV", names=c("Placebo","Vaccine"))
    myboxplot(postvax.B.Brisbane~vaccine, kid, log="y", ylab="HAI titer", main="Cowling", names=c("Placebo","Vaccine"))
mydev.off(file="input/boxplots_rsv_cowling", ext="png")


# make a plot
myfigure (mfrow=c(1,2))
    myboxplot(EIA.log10d14~trt, dat.wide, log="", ylab="",   main="EIA, D14", names=c("Placebo","Vaccine"))
    myboxplot(EIA.log10d14overd0~trt, dat.wide, log="", ylab="",   main="EIA, Fold Change", names=c("Placebo","Vaccine"))
mydev.off(file="input/boxplots_eia", ext="png")


# make a plot
myfigure (mfrow=c(1,2))
    myboxplot(PCA.log10d14~trt, dat.wide, log="", ylab="",   main="PCA, D14", names=c("Placebo","Vaccine"))
    myboxplot(PCA.log10d14overd0~trt, dat.wide, log="", ylab="",   main="PCA, Fold Change", names=c("Placebo","Vaccine"))
mydev.off(file="input/boxplots_pca", ext="png")


# validation using Cowling et al

table(kid$postvax.B.Brisbane)

fits=list()
fits[[1]]=iorw(formula.effect=flu~vaccine+age, formula.mediators=~log2(postvax.B.Brisbane/5), kid, family=binomial(), nboot=0, numCores=1, save.steps=T)
dat.tmp=kid
for (i in 1:9) {
    dat.tmp[dat.tmp$vaccine==0, "postvax.B.Brisbane"]=dat.tmp[dat.tmp$vaccine==0, "postvax.B.Brisbane"]/4
    fits[[1+i]]=iorw(formula.effect=flu~vaccine+age, formula.mediators=~log2(postvax.B.Brisbane/5), dat.tmp, family=binomial(), nboot=0, numCores=1)
}
res=sapply(fits, function(fit) fit$est)
res

myboxplot(postvax.B.Brisbane~vaccine, dat.tmp, log="y")

dat.tmp.2=rbind(dat.tmp, dat.tmp, dat.tmp)
myboxplot(postvax.B.Brisbane~vaccine, dat.tmp.2, log="y")
iorw(formula.effect=flu~vaccine+age, formula.mediators=~log2(postvax.B.Brisbane/5), dat.tmp.2, family=binomial(), nboot=0, numCores=1, save.steps=T)


dat.tmp.3=dat.tmp.2
dat.tmp.3[4,"postvax.B.Brisbane"]=2e-5
myboxplot(postvax.B.Brisbane~vaccine, dat.tmp.3, log="y")
fit.4=iorw(formula.effect=flu~vaccine+age, formula.mediators=~log2(postvax.B.Brisbane/5), dat.tmp.3, family=binomial(), nboot=0, numCores=1, save.steps=F)


dat.tmp.4=kid
myboxplot(postvax.B.Brisbane~vaccine, dat.tmp.4, log="y")
fit.5=iorw(formula.effect=flu~vaccine+age, formula.mediators=~log2(postvax.B.Brisbane/5), dat.tmp.4, family=binomial(), nboot=0, numCores=1, save.steps=T, verbose=F)
print(fit.5)



#### looking at fold change as mediator
dat.wide$tmp=ifelse(dat.wide$EIA.log10d14overd0>1, dat.wide$EIA.log10d14overd0, 0)
myboxplot(tmp~trt, dat.wide, log="")
fit.1=iorw(y2~trt+riskScore.mat.endpoint2+vacc2birth, ~tmp, dat.wide, family=binomial(), nboot=0, numCores=1, save.steps=T); fit.1
#fit.1=iorw(y2~trt+riskScore.mat.endpoint2+vacc2birth, ~tmp, dat.wide, family=binomial(), nboot=1000, numCores=20, save.steps=T); fit.1
##### boot.perc:
#          total        ve    direct  indirect        prop
#2.5%  0.1230524 0.5169078 0.1808655 0.3447321 -0.05110825
#97.5% 0.4830922 0.8769476 1.0134565 1.0665909  1.01252636


dat.wide$tmp=ifelse(dat.wide$EIA.log10d14overd0>log10(4), dat.wide$EIA.log10d14overd0-log10(4), 0)
myboxplot(tmp~trt, dat.wide, log="")
myboxplot(EIA.log10d14overd0~trt, dat.wide, log="")
fit.2=iorw(y2~trt+riskScore.mat.endpoint2+vacc2birth, ~tmp, dat.wide, family=binomial(), nboot=0, numCores=1, save.steps=T); fit.2


#### looking at fold change as mediator
dat.wide$tmp=ifelse(dat.wide$EIA.log10d14>3.75, dat.wide$EIA.log10d14, 0)
myboxplot(tmp~trt, dat.wide, log="")
fit.1=iorw(y2~trt+riskScore.mat.endpoint2+vacc2birth, ~tmp, dat.wide, family=binomial(), nboot=0, numCores=1, save.steps=T); fit.1
