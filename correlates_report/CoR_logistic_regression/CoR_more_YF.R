# @ by Youyi Fong
# start R in this folder or make sure working directory is here
rm(list=ls())    
save.results.to="input"    
library(kyotil)
library(osDesign)
library(aucm)
library(RSVcorr)
library(survey)
library(Hmisc)# wtd.mean
options(survey.lonely.psu = "adjust")# options for svyglm to deal with strata with counts of 1
        
    
summary(subset(dat.wide,trt==1)$EIA.log10d14overd0)
mean(subset(dat.wide,trt==1)$EIA.log10d14overd0>log10(39.1), na.rm=T)

with(dat.wide, table(y1, y2, trt))

with(dat.wide, table(y1, !is.na(EIA.log10d14overd0), trt))
with(dat.wide, table(y2, !is.na(EIA.log10d14overd0), trt))
with(dat.wide, table(y3, !is.na(EIA.log10d14overd0), trt))

with(dat.wide, table(ppefmfl, trt))

with(dat.wide, table(!is.na(EIA.log10d14overd0), trt))


###################################################################################################
# for revision


####################################
# fit a linear regression model of infant age at disease outcome vs. antibody level, in endpoint cases, separately by vaccine vs. placebo. 
# And estimating the treatment arm x antibody level interaction test would address the question of whether vaccination delayed RSV infection. 
# this table is not incorporated into cor report

tab=NULL
for (i in 1:2) {
    tab=rbind(tab, mysapply (assays,function(a) {
        myprint(a, i)
        fit.1=lm(as.formula(paste0("p.1st.epi.st.age ~ ", a, ".log10d14overd0")), dat.wide[dat.wide[["y"%.%i]]==1 & dat.wide$trt==1,])
        fit.2=lm(as.formula(paste0("p.1st.epi.st.age ~ ", a, ".log10d14overd0")), dat.wide[dat.wide[["y"%.%i]]==1 & dat.wide$trt==0,])
        fit.3=lm(as.formula(paste0("p.1st.epi.st.age ~ ", a, ".log10d14overd0*trt")), dat.wide[dat.wide[["y"%.%i]]==1,])
        tab=getFormattedSummary(list(fit.1,fit.2,fit.3), robust=F, type=6)[-1,]
        c(tab[1,1], tab[1,2], tab[3,3])
    }))
}

colnames(tab)=c("vaccine","placebo","itxn")
mytex(tab, file=paste0("input/rev1"), align="c",    
    add.to.row=list(list(0,4), 
        c("       \n \\multicolumn{4}{l}{Endpoint 1} \\\\ \n",
          "\\hline\n \\multicolumn{4}{l}{Endpoint 2}\\\\ \n"
         )
    ))



dat.wide.vacc=subset(dat.wide, trt==1)
dat.wide.plac=subset(dat.wide, trt==0)

table(subset(dat.wide, trt==0)$riskScore.mat.endpoint1)
table(subset(dat.wide, trt==1)$riskScore.mat.endpoint1)

myfigure(mfrow=c(2,2))
    for (a in assays) myboxplot(as.formula(a%.%".log10d14overd0~riskScore.mat.endpoint1"), dat.wide.vacc, test="w", main=a)
mydev.off(file="input/boxplot_foldrise_riskscore")

myfigure(mfrow=c(2,2))
    for (a in assays) myboxplot(as.formula(a%.%".log10d14~riskScore.mat.endpoint1"), dat.wide.vacc, test="w", main=a)
mydev.off(file="input/boxplot_d14_riskscore")




# XX randomized maternal participants did not qualify for the per-protocol efficacy cohort. 

with(dat.wide, table(ppefmfl))

# Supplementary Table xx summarizes the distribution of the 11 maternal enrollment factors used for developing a maternal enrollment risk score, comparing maternal participants who did vs. who did not qualify for the per-protocol efficacy cohort."


###################################################################################################
# Table 1 for the manuscript

digit=0
markers=c((outer(assays[c(3,4,1,2)],".log10"%.%times[c(1,2,3)], paste0))); markers; names(markers)=markers
gmean = function(x, x.ctrl=NULL, na.rm=TRUE, digit=0, wt=rep(1,length(x)), wt.ctrl=rep(1,length(x.ctrl))) {
    n=sum(!is.na(x))
    if(is.null(x.ctrl)) {
        mean=wtd.mean(x, na.rm=na.rm, weights=wt) 
        sd=wtd.var(x, na.rm=na.rm, weights=wt)^.5
    } else {
        # contrast with non-cases
        mean=wtd.mean(x, na.rm=na.rm, weights=wt) - wtd.mean(x.ctrl, na.rm=na.rm, weights=wt.ctrl) 
        sd=  (wtd.var(x, na.rm=na.rm, weights=wt) + wtd.var(x.ctrl, na.rm=na.rm, weights=wt.ctrl))^.5
    }
    tmp = 10^(c(mean, mean-1.96*sd, mean+1.96*sd))
    paste0(formatDouble(tmp[1], digit, remove.leading0=F), " (", formatDouble(tmp[2], digit, remove.leading0=F), "-", formatDouble(tmp[3], digit, remove.leading0=F), ")")
}


for (.trt in c(1,0)) {
#.trt=1; a=markers[1]; i=1
    tab=mysapply (markers, function(a) {
        if (startsWith(a, "RSV")) {
            dat=dat.wide.v
        } else {
            dat=dat.wide
            dat$wt=1 # for wtd.mean
        }
        at.b=endsWith(a, "d0")
        
        out=c()
        
        # non-case
        dat.tmp=subset(dat, trt==.trt & y1==0)    
        m.ctrl=dat.tmp[[a]]
        # m.1.ctrl have baseline subtracted
        if(!at.b) m.1.ctrl=m.ctrl-dat.tmp[[sub("cord","d0",sub("d14","d0",a))]] 
        wt.ctrl=dat.tmp$wt
        out=c(sum(!is.na(m.ctrl)),  gmean(m.ctrl,wt=dat.tmp$wt),  if(!at.b) gmean(m.1.ctrl,wt=dat.tmp$wt) else NA)
        
        # y1 and y2
        for (i in 1:2) {
            dat.tmp=subset(dat, trt==.trt & if (i==1) y1==1 else y2==1)
            m=dat.tmp[[a]]
            # m.1 have baseline subtracted
            if(!at.b) m.1=m-dat.tmp[[sub("cord","d0",sub("d14","d0",a))]]
            
#            # redefine controls if needed
#            if (i==1) {
#                dat.tmp=subset(dat, trt==.trt & y1==0)    
#                m.ctrl=dat.tmp[[a]]
#                if(!at.b) m.1.ctrl=m.ctrl-dat.tmp[[sub("cord","d0",sub("d14","d0",a))]] # subtract baseline
#            } else {
#                dat.tmp=subset(dat, trt==.trt & y2==0)    
#                m.ctrl=dat.tmp[[a]]
#                if(!at.b) m.1.ctrl=m.ctrl-dat.tmp[[sub("cord","d0",sub("d14","d0",a))]] # subtract baseline
#            }
            
            out=c(out, 
                  sum(!is.na(m)),  gmean(m,wt=dat.tmp$wt),                                  if(!at.b) gmean(m.1,wt=dat.tmp$wt) else NA,
                                   gmean(m, m.ctrl,digit=1,wt=dat.tmp$wt,wt.ctrl=wt.ctrl),  if(!at.b) gmean(m.1, m.1.ctrl,digit=1,wt=dat.tmp$wt,wt.ctrl=wt.ctrl) else NA)# substract non-case
                                   
        }    
        
        out    
    })
    print(tab)
    mywrite.csv(tab, file=paste0("input/tab1_", ifelse(.trt==1,"vacc","plac")), row.names=T)   
}       

#with(dat.wide, table(y1, y2, trt))



##########################
# gam and segmented models

myfigure(mfrow=c(1,2))
    fit.1 <- mgcv::gam(y1~s(EIA.log10d14overd0),data=subset(dat.wide,trt==1), family=binomial); plot(fit.1, main="y1")
    fit <-   mgcv::gam(y2~s(EIA.log10d14overd0),data=subset(dat.wide,trt==1), family=binomial); plot(fit, main="y2")
mydev.off(file=paste0("input/EIA.log10d14overd0_gam"), ext="png")

## this takes some time
#fit = chngptm(y2~riskScore.mat.endpoint2+vacc2birth, ~EIA.log10d14overd0, data=subset(dat.wide,trt==1), family="binomial", type="hinge", search.bound=1000, verbose=1, var.type="bootstrap")
#save(fit, file="input/hinge.Rdata")

load(file="input/hinge.Rdata")

fit.step = chngptm(y2~riskScore.mat.endpoint2+vacc2birth, ~EIA.log10d14overd0, data=subset(dat.wide,trt==1), family="binomial", type="step", search.bound=1e5, verbose=1, var.type="none", est.method="grid")

myfigure(mfrow=c(1,2))
    plot(fit, which=1, main="Hinge Model")
    plot(fit.step, which=1, main="Step Model")
mydev.off(file=paste0("input/EIA.log10d14overd0_thresholdmodel"), ext="png")


#########################################################################################
# adjust baseline

myfigure(mfrow=c(1,2))
    plot(EIA.log10d14overd0~EIA.d0, subset(dat.wide, trt==1), log="x", main="Vaccine")
    plot(EIA.log10d14overd0~EIA.d0, subset(dat.wide, trt==0), log="x", main="Placebo")
mydev.off(file=paste0(save.results.to, "/scatterplot_EIA_foldchange_d0"), ext="png")    


# when 
fit=glm(formula=y2~log10(EIA.d0), subset(dat.wide, trt==1), family="binomial"); summary(fit)
fit=glm(formula=y2~EIA.log10d14overd0, subset(dat.wide, trt==1), family="binomial"); summary(fit)
fit=glm(formula=y2~EIA.log10d14overd0+log10(EIA.d0), subset(dat.wide, trt==1), family="binomial"); summary(fit)



#########################################################################################
# boxplots by cases and controls

# Obj 1 boxplot for fold change

for (i in 1:3) {
for (.trt in 1:0) {    
    myfigure(mfcol=c(4,6))    
        for (a in assays) {
            dat = subset(dat.wide.v, trt==.trt)
            myboxplot(as.formula(a%.%".d0~y"%.%i), dat, main=paste0(a), ylab="D0", names=c("controls","cases"), log="y")
            myboxplot(as.formula(a%.%".d14~y"%.%i), dat, main=paste0(a), ylab="D14", names=c("controls","cases"), log="y")
            myboxplot(as.formula(a%.%".log10d14overd0~y"%.%i), dat, main=paste0(a), ylab="log10 fold change", names=c("controls","cases"))
            myboxplot(as.formula(a%.%".cord~y"%.%i), dat, main=paste0(a), ylab="cord blood", names=c("controls","cases"), log="y")
        }    
    
        for (a in assays[3:4]) {
            dat = subset(dat.wide, trt==.trt)            
            myboxplot(as.formula(a%.%".d0~y"%.%i), dat, main=paste0(a), ylab="D0", names=c("controls","cases"), log="y")
            myboxplot(as.formula(a%.%".d14~y"%.%i), dat, main=paste0(a), ylab="D14", names=c("controls","cases"), log="y")
            myboxplot(as.formula(a%.%".log10d14overd0~y"%.%i), dat, main=paste0(a), ylab="log10 fold change", names=c("controls","cases"))
            myboxplot(as.formula(a%.%".cord~y"%.%i), dat, main=paste0(a), ylab="cord blood", names=c("controls","cases"), log="y")
        }    
    mydev.off(file=paste0(save.results.to, "/boxplot_by_y",i,"_", ifelse(.trt==1,"vacc","plac")), ext="png")    
}
}

for (t in c("log10d14overd0", "log10d14", "log10d0", "log10cord")){
    for (i in 1:3) {
    for (.trt in 0:1) {    
        myfigure(mfcol=c(1,4))    
            for (a in assays[1:2]) {
                dat = subset(dat.wide.v, trt==.trt)
                myboxplot(as.formula(a%.%"."%.%t%.%"~y"%.%i), dat, main=paste0(a), ylab="log10 fold change", names=c("controls","cases"))
            }        
            for (a in assays[3:4]) {
                dat = subset(dat.wide, trt==.trt)            
                myboxplot(as.formula(a%.%"."%.%t%.%"~y"%.%i), dat, main=paste0(a), ylab="log10 fold change", names=c("controls","cases"))
            }    
        mydev.off(file=paste0(save.results.to, "/boxplot_",t,"_by_y",i,"_", ifelse(.trt==1,"vacc","plac")), ext="png")    
    }
    }
}


# interaction plots
for (i in 1:3) {
for (.trt in 1:1) {    
# i=1; .trt=1
    myfigure(mfcol=c(1,4))    
        for (a in assays[1:4]) {
            dat = subset(dat.wide, trt==.trt)[c(a%.%".log10d0",a%.%".log10d14")]
            dat=dat[complete.cases(dat),]
            my.interaction.plot(dat, x.ori = 0, xaxislabels = rep("", 2), cex.axis = 1, add = FALSE, xlab = "", ylab = "log10 titer", pcol = NULL, lcol = 1, main=a)
            axis(1, 1:2, c("D0","D14"))
            dat = subset(dat.wide, trt==.trt & dat.wide[["y"%.%i]]==1)[c(a%.%".log10d0",a%.%".log10d14")]
            dat=dat[complete.cases(dat),]
            my.interaction.plot(dat, x.ori = 0, xaxislabels = rep("", 2), cex.axis = 1, add = T, pcol = NULL, lcol = 2)
        }    
    mydev.off(file=paste0(save.results.to, "/interactionplot_d0_d14_by_y",i,"_", ifelse(.trt==1,"vacc","plac")), ext="png")    
}
}



# cutoff based on looking at the boxplots
dat.wide$EIA.d0pos=ifelse(dat.wide$EIA.log10d0>2.5,1,0)
dat.wide$PCA.d0pos=ifelse(dat.wide$PCA.log10d0>1,1,0)
with(subset(dat.wide, y2==1 & trt==1), table(EIA.d0pos))
with(subset(dat.wide, y2==1 & trt==1), table(PCA.d0pos))
with(subset(dat.wide, y1==1 & trt==1), table(EIA.d0pos))
with(subset(dat.wide, y1==1 & trt==1), table(PCA.d0pos))

# separate CoR analyses for baselinen pos and neg
for (y in c("y1", "y2")) {
for (x in c("d14overd0","d14")) {
for (m in c("EIA","PCA")) {
    
    formula=paste0(y, "~", m, ".log10", x)
    # add maternal risk scores
    if (y=="y1") formula=paste0(formula,"+riskScore.mat.endpoint1") else formula=paste0(formula,"+riskScore.mat.endpoint2")
    # adj vacc2birth for maternal timepoints
    if(!contain(x,"cord")) formula=paste0(formula,"+vacc2birth")
    formula=as.formula(formula); print(formula)
    
    # baseline pos+neg
    fit.1=glm(formula=formula, subset(dat.wide, trt==1), family="binomial")
    # baseline pos
    fit.2=glm(formula=formula, subset(dat.wide, trt==1 & dat.wide[[m%.%".d0pos"]]==1), family="binomial")
    # baseline neg
    fit.3=glm(formula=formula, subset(dat.wide, trt==1 & dat.wide[[m%.%".d0pos"]]==0), family="binomial")
    # interaction
    fit.4=glm(formula=update(formula, as.formula(paste0("~.+", m, ".log10", x, "*", m, ".d0pos"))), subset(dat.wide, trt==1), family="binomial")
    # main
    fit.5=glm(formula=update(formula, as.formula(paste0("~.+", m, ".d0pos"))), subset(dat.wide, trt==1), family="binomial")
    
    tab=getFormattedSummary(list(fit.1, fit.2, fit.3, fit.5, fit.4), robust=F, exp=T, type=6)
    colnames(tab)=c("D0 pos+neg","D0 pos","D0 neg","main","itxn")
    tab=tab[c(3,4,2,5,6),]; tab
    mytex(tab, file=paste0("tables/CoR_", y, "_trt",.trt, "_",m, "_",x, "_baselinepos"), align="c", input.folder=save.results.to)
        
}
}
}




#########################################################################################
# why EIA and PCA together so significant

dat.tmp=subset(dat.wide, trt==1)
dat.tmp.2=subset(dat.wide.v, trt==1)

# scatterplot of EIA and PCA
myfigure(mfrow=c(1,2))
    plot(EIA.d14~PCA.d14, dat.tmp, log="xy", col=dat.tmp$y1+1, pch=ifelse(dat.tmp$y1==1, 2, 1), lwd=ifelse(dat.tmp$y1==1, 2, 1), main="Phase 1")
    plot(EIA.d14~PCA.d14, dat.tmp.2, log="xy", col=dat.tmp.2$y1+1, pch=ifelse(dat.tmp.2$y1==1, 2, 1), lwd=ifelse(dat.tmp.2$y1==1, 2, 1), main="Phase 2")
mydev.off(file=paste0(save.results.to, "/scatterplot_eia_pca_d14_vacc"), ext="png")    


dat.tmp$tmp=with(dat.tmp, log10(EIA.d14)-log10(PCA.d14))
fit=glm(formula=y1~tmp+riskScore.mat.endpoint1+vacc2birth, dat.tmp, family=binomial); summary(fit)

# glm fit to phase 1 data
dat.tmp=subset(dat.wide, trt==1 & !is.na(EIA.d14) & !is.na(PCA.d14))
fit.0=glm(formula=y1~riskScore.mat.endpoint1+vacc2birth, dat.tmp, family=binomial); summary(fit.0)
fit.1=glm(formula=y1~riskScore.mat.endpoint1+vacc2birth+log10(EIA.d14), dat.tmp, family=binomial); summary(fit.1)
fit.2=glm(formula=y1~riskScore.mat.endpoint1+vacc2birth+log10(PCA.d14), dat.tmp, family=binomial); summary(fit.2)
fit.3=glm(formula=y1~riskScore.mat.endpoint1+vacc2birth+log10(EIA.d14)+log10(PCA.d14), dat.tmp, family=binomial); summary(fit.3)
anova(fit.0, fit.1, test="Chisq")
anova(fit.0, fit.2, test="Chisq")
anova(fit.0, fit.3, test="Chisq")

var.ind=4:5
stat=coef(fit.3)[var.ind] %*% solve(vcov(fit.3,robust=F)[var.ind,var.ind]) %*% coef(fit.3)[var.ind] # robust=F is needed for tps fit, not used for tps fit
pchisq(stat, 4, lower.tail = FALSE)


fit=tps.rsv(formula=y1~EIA.d14+PCA.d14 +riskScore.mat.endpoint1+vacc2birth, trt.grp=1)
summary(fit)
var.ind=2:3
stat=coef(fit)[var.ind] %*% solve(vcov(fit,robust=F)[var.ind,var.ind]) %*% coef(fit)[var.ind] # robust=F is needed for tps fit, not used for tps fit
pchisq(stat, 4, lower.tail = FALSE)

fit=tps.rsv(formula=y1~EIA.d14 +riskScore.mat.endpoint1+vacc2birth, trt.grp=1); summary(fit)
fit=tps.rsv(formula=y1~PCA.d14 +riskScore.mat.endpoint1+vacc2birth, trt.grp=1); summary(fit)
fit=tps.rsv(formula=y1~EIA.d14+PCA.d14 +riskScore.mat.endpoint1+vacc2birth, trt.grp=1); summary(fit)



dstrat<-svydesign(id=~1,strata=~stratum, weights=~wt, data=subset(dat.wide.v, trt==.trt))
fit=svyglm(formula=y1~EIA.d14+PCA.d14 +riskScore.mat.endpoint1+vacc2birth, design=dstrat, family="binomial")
summary(fit)
var.ind=2:3
stat=coef(fit)[var.ind] %*% solve(vcov(fit,robust=F)[var.ind,var.ind]) %*% coef(fit)[var.ind] # robust=F is needed for tps fit, not used for tps fit
pchisq(stat, 4, lower.tail = FALSE)



########################################################################################
# fitting EIA/PCA with phase 2 data 

.trt=1

fit.1=tps.rsv(formula=y1~log10(PCA.log10d14overd0)+riskScore.mat.endpoint1+vacc2birth, trt.grp=.trt); summary(fit.1)
fit.2=glm    (formula=y1~log10(EIA.d14)+riskScore.mat.endpoint1+vacc2birth, subset(dat.wide, trt==.trt), family=binomial); summary(fit.2)

dstrat<-svydesign(id=~1,strata=~stratum, weights=~wt, data=subset(dat.wide.v, trt==.trt))
fit.3=svyglm(y1~log10(PCA.log10d14overd0)+riskScore.mat.endpoint1+vacc2birth, design=dstrat, family="binomial"); summary(fit.3)



########################################################################################
#

library(aucm)
library(RSVcorr)
library(WeightedROC)

with(subset(dat.wide, trt==1), fast.auc(riskScore.mat.endpoint2, y2))
# 0.6027337
with(subset(dat.wide.v, trt==1), fast.auc(riskScore.mat.endpoint2, y2)) # the weighted version is .553 or 0.5886225?
# 0.5810772
with(subset(dat.wide.v, trt==1), WeightedAUC(WeightedROC(riskScore.mat.endpoint2, y2, weight = wt)))
# 0.5886225

with(subset(dat.wide, trt==1), fast.auc(child5, y2))
# 0.6451847
with(subset(dat.wide.v, trt==1), fast.auc(child5, y2)) 
# 0.6232304
with(subset(dat.wide.v, trt==1 & !is.na(child5)), WeightedAUC(WeightedROC(as.integer(child5=="Y"), y2, weight = wt)))
# 0.6393662


########################################################################################
#


with(subset(dat.wide, y1==1), table(itteffl=="Y", ppefmfl=="Y"))


with(subset(dat.wide, y1==1), table(ppefifl=="Y", ppefmfl=="Y"))



with(dat.wide, table(y1, ppefmfl=="Y"))
with(dat.wide, table(y1, trt, ppefifl=="Y"))
with(dat.wide, table(trt, ppefifl=="Y", y1))

with(dat.wide, table(y2, ppefmfl=="Y"))
with(dat.wide, table(y2, ppefifl=="Y"))
with(dat.wide, table(trt, ppefifl=="Y", y2))



with(dat.wide, table(ppimmfl=="Y", ppefifl=="Y"))

with(dat.wide, table(trt, y1, ppefifl=="Y"))
with(dat.wide.v, table(trt, y1, ppefifl=="Y"))

with(dat.wide, table(trt, y2, ppefifl=="Y"))
with(dat.wide.v, table(trt, y2, ppefifl=="Y"))


with(subset(dat.wide, trt==1), sum(y2))

with(subset(dat.wide, trt==0), sum(y1))
with(subset(dat.wide, trt==0), sum(y2))

with(dat.wide, sum(y1))
with(dat.wide, sum(y2))


table(dat.wide$trt, dat.wide$y1)
table(dat.wide.v$trt, dat.wide.v$y1)



summary(subset(dat.wide, ppimmfl=="Y"))
summary(subset(dat.wide, ppimifl=="Y"))


with(dat.wide, table(ppimmfl=="Y", ppimifl=="Y", y1))
with(dat.wide, table(ppimmfl=="Y", ppimifl=="Y", y2))


with(dat.wide, table(is.na(EIA.log10d14overd0), y2))
with(dat.wide, table(is.na(EIA.cord), y2))



with(dat.wide, table(ppimmfl=="Y", ittimfl=="Y"))


with(subset(dat.wide, y1=="1"), table(ppimmfl=="Y", ppimifl=="Y"))
with(subset(dat.wide, y1=="1"), table(ppimmfl=="Y", ppimifl=="Y", !is.na(EIA.log10d14overd0)))
with(subset(dat.wide, y1=="1"), table(ppimifl=="Y", !is.na(EIA.log10d14overd0), ppimmfl=="Y" ))

with(subset(dat.wide, y1=="1"), table(!is.na(EIA.log10d14overd0), ppimmfl=="Y", trt ))
with(subset(dat.wide, y2=="1"), table(!is.na(EIA.log10d14overd0), ppimmfl=="Y", trt ))

with(subset(dat.wide, y1=="1"), table(!is.na(EIA.log10d14overd0), ppimifl=="Y", trt ))
with(subset(dat.wide, y2=="1"), table(!is.na(EIA.log10d14overd0), ppimifl=="Y", trt ))

with(subset(dat.wide, y1=="1"), table(!is.na(EIA.log10d14overd0), ppimifl=="Y", trt ))
with(subset(dat.wide, y2=="1"), table(!is.na(EIA.log10d14overd0), ppimifl=="Y" & ppimmfl=="Y", trt ))


with(subset(dat.wide, y1=="1"), table(!is.na(EIA.log10d14overd0)))


with(subset(dat.wide, y1=="1"), table(ppimmfl=="Y", ppimifl=="Y"))



with(subset(dat.wide, y1=="1"), table(ppefmfl=="Y", ppimmfl=="Y"))

with(subset(dat.wide, y1=="1"), table(ppefifl=="Y", ppimifl=="Y"))



##
adefepip$trt=dat.wide$trt[match(adefepip$SCRMSUBJID, dat.wide$pair.id%.%"-P")]

nrow(subset(adefepip, (endpt3_1 == "Y" & endpt3_2!= "Y" &  endpt3_3!= "Y" &  epi_st_age <= 90 & trt==1))) #14
nrow(subset(adefepip, (endpt3_1 != "Y" & endpt3_2== "Y" &  endpt3_3!= "Y" &  epi_st_age <= 90 & trt==1))) #0
nrow(subset(adefepip, (endpt3_1 != "Y" & endpt3_2!= "Y" &  endpt3_3== "Y" &  epi_st_age <= 90 & trt==1))) #14

nrow(subset(adefepip, (endpt3_1 == "Y" & endpt3_2== "Y" &  endpt3_3== "Y" &  epi_st_age <= 90 & trt==1)))#10
nrow(subset(adefepip, (endpt3_1 == "Y" & endpt3_2== "Y" &  endpt3_3!= "Y" &  epi_st_age <= 90 & trt==1)))#4
nrow(subset(adefepip, (endpt3_1 == "Y" & endpt3_2!= "Y" &  endpt3_3== "Y" &  epi_st_age <= 90 & trt==1)))#10
nrow(subset(adefepip, (endpt3_1 != "Y" & endpt3_2== "Y" &  endpt3_3== "Y" &  epi_st_age <= 90 & trt==1)))#0


nrow(subset(adefepip, (endpt3_1 == "Y" & endpt3_2!= "Y" &  endpt3_3!= "Y" &  epi_st_age <= 90 & trt==0))) #3
nrow(subset(adefepip, (endpt3_1 != "Y" & endpt3_2== "Y" &  endpt3_3!= "Y" &  epi_st_age <= 90 & trt==0))) #0
nrow(subset(adefepip, (endpt3_1 != "Y" & endpt3_2!= "Y" &  endpt3_3== "Y" &  epi_st_age <= 90 & trt==0))) #6

nrow(subset(adefepip, (endpt3_1 == "Y" & endpt3_2== "Y" &  endpt3_3== "Y" &  epi_st_age <= 90 & trt==0)))#23
nrow(subset(adefepip, (endpt3_1 == "Y" & endpt3_2== "Y" &  endpt3_3!= "Y" &  epi_st_age <= 90 & trt==0)))#3
nrow(subset(adefepip, (endpt3_1 == "Y" & endpt3_2!= "Y" &  endpt3_3== "Y" &  epi_st_age <= 90 & trt==0)))#15
nrow(subset(adefepip, (endpt3_1 != "Y" & endpt3_2== "Y" &  endpt3_3== "Y" &  epi_st_age <= 90 & trt==0)))#1


with(subset(dat.wide, trt==0), table(y1, y2))


# 
subset(dat.wide, EIA.log10d14>4 & y2==1, select=c(EIA.log10d0, EIA.log10d14, EIA.log10d14overd0))

subset(dat.wide, EIA.log10d0>3.5 & y2==1, select=c(EIA.log10d0, EIA.log10d14, EIA.log10d14overd0))



with(subset(dat.wide, trt==1), fast.auc(EIA.log10d14overd0, y2))
with(subset(dat.wide, trt==1), fast.auc(PCA.log10d14overd0, y2))
with(subset(dat.wide, trt==1), fast.auc(RSVA.log10d14overd0, y2))
with(subset(dat.wide, trt==1), fast.auc(RSVB.log10d14overd0, y2))

library(WeightedROC)
with(subset(dat.wide.v, trt==1), WeightedAUC(WeightedROC(1-RSVA.log10d14overd0, y2, wt)))
with(subset(dat.wide.v, trt==1), WeightedAUC(WeightedROC(1-RSVB.log10d14overd0, y2, wt)))
