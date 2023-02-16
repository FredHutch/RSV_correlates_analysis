rm(list=ls())
# @ by Youyi Fong
# paths that may need to be changed to reproduce results
if (file.exists("c:/_checkReproducibility")) setwd("D:/gdrive/RSV/R")
save.results.to="D:/gdrive/RSV/R/RSVcorrelatesAnalysis/correlates_report/input"
    
library(RSVcorr) 
library(kyotil)
library(osDesign) # tps
library(aucm) # fastauc
library(survey) # svyglm
options(survey.lonely.psu = "adjust")# options for svyglm to deal with strata with counts of 1

library(doParallel)
numCores=30 # volta

dat.wide.vacc=subset(dat.wide, trt==1)
f1=~riskScore.mat.endpoint2 + vacc2birth
f2=~riskScore.mat.endpoint2


########################################################################################
# 

times=c("log10d0","log10d14","log10d14overd0","log10cord")
markers=c(t(outer(assays,"."%.%times, paste0))); markers; names(markers)=markers    

seeds=1:1e4
out=mclapply(seeds, mc.cores = numCores, FUN=function(seed) {   
    set.seed(seed)    
    new.idx=sample(1:nrow(dat.wide.vacc))
    pvals=sapply (markers, function(a) {
        dat.wide.vacc[[a]]=dat.wide.vacc[[a]][new.idx]# permutate marker only
        f=if(!endsWith(a, "cord")) f1 else f2
        fit=glm(update(f, as.formula("y2~.+"%.%a)), dat.wide.vacc, family = binomial())
        last(c(summary(fit)$coef))
    })        
}) 
out=do.call(cbind, out)


ref.dist=apply(out, 2, min)
mean(ref.dist<0.00627)
# xx if d0, d14, fold change, cord are adjusted
# smaller if d0, d14, fold change are adjusted
# even smaller if d14, fold change are adjusted


# The table CoR_obj1a_y2_trt1_padj.tex is probably made manually partially based on results from this file
