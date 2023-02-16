# start R in this folder or make sure working directory is here
rm(list=ls())    
save.results.to="input"    
library(kyotil)
library(osDesign)
library(survey)
library(RSVcorr)
# options for svyglm to deal with strata with counts of 1
options(survey.lonely.psu = "adjust")


#########################################################################################
# CoR Object 2a
# putting different assays in the same models, one model for each time point

var.ind=2:5
for (y in c("y1","y2","y3")) {
for (.trt in c(1,0)) {
for(ind in 1:2) { # 1: not controlling for baseline, 2: controlling for baseline RSVA/RSVB average
#y="y1"; .trt=1; ind=2; x="log10d14overd0"
#myprint(y, .trt, ind)
    
    tab=NULL
    fits=list()
    for (x in times[c(1,2,4,3)]) {
        markers=c(t(outer(assays,"."%.%x, paste0)))    
        markers=if(contain(x,"log10")) markers else "log10("%.%markers%.%")"
        formula = paste0(y, "~", concatList(markers,"+"))
        # add maternal risk scores
        if (y=="y1") formula=paste0(formula,"+riskScore.mat.endpoint1") else formula=paste0(formula,"+riskScore.mat.endpoint2")
        # adj vacc2birth for maternal timepoints
        if(!contain(x,"cord")) formula=paste0(formula,"+vacc2birth")
    
        if (ind==2 & x!="d0") {
            dat.wide.v$log10.RSVAB.d0 = log10((dat.wide.v$RSVA.d0+dat.wide.v$RSVB.d0)/2)
            formula=paste0(formula, "+log10.RSVAB.d0")
        } 
        formula=as.formula(formula); formula
        
        if (y=="y1") {
            # tps regression
            fit=tps.rsv(formula=formula, trt.grp=.trt)
        } else {
            # ipw
            dstrat<-svydesign(id=~1,strata=~stratum, weights=~wt, data=subset(dat.wide.v, trt==.trt))
            fit=svyglm(formula, design=dstrat, family="binomial")
        }                                          
        fits[[x]]=fit
    }        
    
    if (ind==1) {
        tab=getFormattedSummary(fits, robust=F, exp=T, rows=var.ind)# rows is important for being able to extract info from both cord and d0 etc
        colnames(tab)= ifelse(contain(colnames(tab),"log10"),colnames(tab),"log10("%.%colnames(tab)%.%")")
        rownames(tab)= sub(".d0)","",rownames(tab))
        rownames(tab)= sub("log10\\(","",rownames(tab))
        tab    
        mytex(tab, file=paste0("tables/CoR_obj2a",ind,"_", y, "_trt",.trt), align="c", input.folder=save.results.to)
    
    } else {
        tab=sapply(1:length(fits), simplify="array", function(i) {
            out=getFormattedSummary(fits[i], robust=F, exp=T)
            out=out[contain(rownames(out),"log10"),]
            if(i==1) c(out,"RSVA/B.log10d0"=NA) else out
        })
        colnames(tab)= ifelse(contain(times[c(1,2,4,3)],"log10"),times[c(1,2,4,3)],"log10("%.%times[c(1,2,4,3)]%.%")")
        rownames(tab)= sub(".d0)","",rownames(tab))
        rownames(tab)= sub("log10\\(","",rownames(tab))
        tab    
        mytex(tab, file=paste0("tables/CoR_obj2a",ind,"_", y, "_trt",.trt), align="c", input.folder=save.results.to)    
    }
    
        
    # multitesting adjustment
    # get Wald test p value
    pvals=sapply(fits, function(fit) {
        stat=coef(fit)[var.ind] %*% solve(vcov(fit,robust=F)[var.ind,var.ind]) %*% coef(fit)[var.ind] # robust=F is needed for tps fit, not used for tps fit
        pchisq(stat, 4, lower.tail = FALSE)
    }) 
    tab.2=cbind(pvals, Holm=p.adjust(pvals, method="holm"), BH=p.adjust(pvals, method="BH")); tab.2
    mytex(tab.2, file=paste0("tables/CoR_obj2a",ind,"_", y, "_trt",.trt, "_pvals"), align="c", input.folder=save.results.to)
    
}
}
}


# obj2a3: only EIA and PCA, phase 1 data
var.ind=2:3
for (y in c("y1","y2","y3")) {
for (.trt in c(1,0)) {
#y="y1"; .trt=1; x="d14"
    
    tab=NULL
    fits=list()
    for (x in times[c(1,2,4,3)]) {
        markers=c(t(outer(assays[3:4],"."%.%x, paste0)))    
        markers=if(contain(x,"log10")) markers else "log10("%.%markers%.%")"
        formula = paste0(y, "~", concatList(markers,"+"))
        # add maternal risk scores
        if (y=="y1") formula=paste0(formula,"+riskScore.mat.endpoint1") else formula=paste0(formula,"+riskScore.mat.endpoint2")
        # adj vacc2birth for maternal timepoints
        if(!contain(x,"cord")) formula=paste0(formula,"+vacc2birth")
        fit=glm(formula=formula, subset(dat.wide, trt==.trt), family="binomial")
        fits[[x]]=fit
    }        
    
    tab=getFormattedSummary(fits, robust=F, exp=T, rows=var.ind)# rows is important for being able to extract info from both cord and d0 etc
    colnames(tab)= ifelse(contain(colnames(tab),"log10"),colnames(tab),"log10("%.%colnames(tab)%.%")")
    rownames(tab)= sub(".d0)","",rownames(tab))
    rownames(tab)= sub("log10\\(","",rownames(tab))
    tab    
    mytex(tab, file=paste0("tables/CoR_obj2a3_", y, "_trt",.trt), align="c", input.folder=save.results.to)
    
    # multitesting adjustment
    # get Wald test p value
    pvals=sapply(fits, function(fit) {
        stat=coef(fit)[var.ind] %*% solve(vcov(fit,robust=F)[var.ind,var.ind]) %*% coef(fit)[var.ind] # robust=F is needed for tps fit, not used for tps fit
        pchisq(stat, 4, lower.tail = FALSE)
    }) 
    tab.2=cbind(pvals, Holm=p.adjust(pvals, method="holm"), BH=p.adjust(pvals, method="BH")); tab.2
    mytex(tab.2, file=paste0("tables/CoR_obj2a3_", y, "_trt",.trt, "_pvals"), align="c", input.folder=save.results.to)

    # likelihood ratio test
    pvals=sapply (times[c(1,2,4,3)], function(x) {
        markers=c(t(outer(assays[3:4],"."%.%x, paste0)))            
        
        dat.tmp=subset(dat.wide, trt==.trt)
        dat.tmp=dat.tmp[apply(dat.tmp[markers], 1, function(x) all(!is.na(x))),]
        
        markers=if(contain(x,"log10")) markers else "log10("%.%markers%.%")"
        formula = paste0(y, "~", concatList(markers,"+"))
        # add maternal risk scores
        if (y=="y1") formula=paste0(formula,"+riskScore.mat.endpoint1") else formula=paste0(formula,"+riskScore.mat.endpoint2")
        # adj vacc2birth for maternal timepoints
        if(!contain(x,"cord")) formula=paste0(formula,"+vacc2birth")
        f0=paste0(y, ifelse (y=="y1", "~riskScore.mat.endpoint1", "~riskScore.mat.endpoint2"), ifelse(!contain(x,"cord"), "+vacc2birth", ""))
        
        fit.1=glm(formula=formula, dat.tmp, family="binomial")
        fit.0=glm(formula=f0     , dat.tmp, family="binomial")
        anova(fit.1, fit.0, test="Chisq")$P[2]
    })
    tab.3=cbind(pvals, Holm=p.adjust(pvals, method="holm"), BH=p.adjust(pvals, method="BH")); tab.3
    mytex(tab.3, file=paste0("tables/CoR_obj2a3_", y, "_trt",.trt, "_pvals_lr"), align="c", input.folder=save.results.to)


    
}
}






#########################################################################################
# CoR Object 2b
# putting different time points in the same models, one model for each assay

var.ind=2:4
for (y in c("y1","y2","y3")) {
for (.trt in c(1,0)) {
for(ind in 1:2) { # 1: d0, d14, cord; 2: d0, d14overd0, cord
#y="y1"; .trt=0; ind=1

    tab=NULL
    fits=list()
    for (x in assays) {
        markers=c(t(outer(x, "."%.%times[if(ind==1) 1:3 else c(1,4,3)], paste0)))    
        markers=ifelse(contain(markers,"log10"), markers, "log10("%.%markers%.%")")
        formula = paste0(y, "~", concatList(markers,"+")) 
        # add maternal risk scores
        if (y=="y1") formula=paste0(formula,"+riskScore.mat.endpoint1") else formula=paste0(formula,"+riskScore.mat.endpoint2")
        # vacc2birth not adjusted because cord measurement is adjusted
        formula=as.formula(formula); formula
        
        if(x %in% c("RSVA","RSVB")) {
            # phase 2 data
            if (y=="y1") {
                fit=tps.rsv(formula=formula, trt.grp=.trt)
            } else {
                dstrat<-svydesign(id=~1,strata=~stratum, weights=~wt, data=subset(dat.wide.v, trt==.trt))
                fit=svyglm(formula, design=dstrat, family="binomial")
            }                                          
        } else {
            # phase 1 data
            fit=glm(formula=formula, subset(dat.wide, trt==.trt), family="binomial")
        }
        fits[[x]]=fit
    }        
    
    tab=getFormattedSummary(fits, robust=F, exp=T)[-1,]
    rownames(tab)= sub("RSVA.","",rownames(tab))
    tab
    
    mytex(tab, file=paste0("tables/CoR_obj2b",ind,"_", y, "_trt",.trt), align="c", input.folder=save.results.to)

    # multitesting adjustment
    # get Wald test p value
    pvals=sapply(fits, function(fit) {
        stat=coef(fit)[var.ind] %*% solve(vcov(fit,robust=F)[var.ind,var.ind]) %*% coef(fit)[var.ind] # robust=F is needed for tps fit, not used for tps fit
        pchisq(stat, 4, lower.tail = FALSE)
    }) 
    tab.2=cbind(pvals, Holm=p.adjust(pvals, method="holm"), BH=p.adjust(pvals, method="BH")); tab.2
    mytex(tab.2, file=paste0("tables/CoR_obj2b",ind,"_", y, "_trt",.trt, "_pvals"), align="c", input.folder=save.results.to)
    
}    
}
}
