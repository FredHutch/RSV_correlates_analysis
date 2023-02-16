# start R in this folder or make sure working directory is here
rm(list=ls())    
rerun.marginal.risk.boot=T # change to F if bootstrap results are saved
    
save.results.to="input"    
library(kyotil); stopifnot(packageVersion("kyotil")>="2021.1.4")
library(osDesign)
library(survey)
# options for svyglm to deal with strata with counts of 1
options(survey.lonely.psu = "adjust")
#library(doParallel) # loading this overwrites RSVcorr::times and there is no way to get that back
numCores=30 # volta
library(Hmisc)
library(forestplot); stopifnot(packageVersion("forestplot")=="1.10") # newer ones have a bug such that when in a loop, empty forest plots are generated
library(marginalizedRisk); stopifnot(packageVersion("marginalizedRisk")>="2020.12.17")
library(RSVcorr); stopifnot(packageVersion("RSVcorr")>="2020.10.31")
dat.wide.vacc=subset(dat.wide, trt==1)
dat.wide.plac=subset(dat.wide, trt==0)





#########################################################################################
# time consuming step, to be run on server

# simplest bootstrap, all samples are available (no phase 2), 
marginal.risk.boot=function(formula, marker.name, data, B, ci.type="quantile", numCores=1) {  
    
    # store the current rng state 
    save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
    if (class(save.seed)=="try-error") {set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) }      
    
    ss=quantile(data[[marker.name]], seq(.05,.95,by=0.01))
    n=row(data)
    
    fit.risk=glm(update(formula, as.formula(paste0("~.+",marker.name))), data, family=binomial)
    fit.s=lm(update(formula, as.formula(paste0(marker.name,"~."))), data) 
    prob=marginal.risk.cont(fit.risk, fit.s, data, ss=ss)
    #prob=marginalized.risk(fit.risk, marker.name, data=data.ph2, ss=ss, weights=data.ph2$wt, t=t, categorical.s=F)    

    
    out=mclapply(1:B, mc.cores = numCores, FUN=function(seed) {   
        set.seed(seed)    
        dat.tmp.mrb=data[sample(n, replace=TRUE),]        
        f1=update(formula, as.formula(paste0("~.+",marker.name)))
        f2=update(formula, as.formula(paste0(marker.name,"~.")))
## this is once thought needed
#        environment(f1) <- list2env(list(data=dat.tmp.mrb)) 
#        environment(f2) <- list2env(list(data=dat.tmp.mrb)) 
        fit.risk=glm(f1, dat.tmp.mrb, family=binomial)
        fit.s=    lm(f2, dat.tmp.mrb) 
        marginal.risk.cont(fit.risk, fit.s, dat.tmp.mrb, ss=ss)
    })
    res=do.call(cbind, out)
    
    # restore rng state 
    assign(".Random.seed", save.seed, .GlobalEnv)    
    
    if (ci.type=="quantile") {
        ci.band=t(apply(res, 1, function(x) quantile(x, c(.025,.975))))
    } else {
        stop("only quantile bootstrap CI supported for now")
    }
    
    list(marker=ss, prob=prob, boot=res, lb=ci.band[,1], ub=ci.band[,2])     
}    

if (rerun.marginal.risk.boot) {    
    # do this a machine with at least numCores cpus
    # this takes a few hours to run
    risks=list()    
    for (i in 1:3) {
        if (i==1) {
            formula=y2~riskScore.mat.endpoint2+vacc2birth
            marker.name="EIA.log10d14overd0"
        } else if(i==2) {
            formula=y2~riskScore.mat.endpoint2+vacc2birth
            marker.name="EIA.log10d14"
        } else if(i==3) {
            formula=y2~riskScore.mat.endpoint2+vacc2birth+EIA.log10d0
            marker.name="EIA.log10d14overd0"
        }
        dat.tmp=dat.wide.vacc[!is.na(dat.wide.vacc[[marker.name]]), ] # make sure there is no missing data
        risks[[i]]= marginal.risk.boot(formula, marker.name, dat.tmp, B=1000, numCores=numCores) 
    }    
    if (!dir.exists("input")) dir.create("input")
    save(risks, file="input/risks.Rdata")
} else {
    load(file="input/risks.Rdata")
}


# 26Oct2020      Erika Rudnicki
# to make CI wider, make width bigger and graphwidth larger
theforestplot <- function(group,point.estimates,lower.bounds,upper.bounds,p.values,table.labels,zero.line=1.0,dashed.line=NA,x.ticks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2),
                          decimal.places = 1,fontsize = 1,width.pdf=7,height.pdf=7,graphwidth="auto",...){
  
  plotdata <- structure(
    list(
      mean  = c(NA, point.estimates),
      lower = c(NA, lower.bounds),
      upper = c(NA, upper.bounds)
    ),
    .Names = c("mean", "lower", "upper"),
    row.names = c(NA,-(length(point.estimates) + 1)),
    class = "data.frame"
  )
  
  # remove.redundancy <- function(x){ifelse(!is.na(dplyr::lag(x)) & x==dplyr::lag(x), NA, x)}
  show.decimals <- function(x){format(round(x, decimal.places), nsmall=decimal.places)}
  
  tabletext <- cbind(
    # c(table.labels[1], remove.redundancy(as.character(cohort))),
      c(table.labels[1], as.character(group))
    , c(paste0(table.labels[2]), paste0(sapply(point.estimates, show.decimals), " (", sapply(lower.bounds, show.decimals), ", ", sapply(upper.bounds, show.decimals), ")"))
    , c("P-value", formatDouble(p.values, 3, remove.leading0 =F))
  )
    
  forestplot(tabletext,plotdata,is.summary = FALSE,col = fpColors(box = "darkblue",line = "darkblue",summary = "royalblue",zero="black"),    
    graph.pos = 4,graphwidth = graphwidth,hrzl_lines = list("2" = gpar(lty=1)),zero = zero.line,lwd.zero = 0.5,lwd.ci = 0.5,lwd.xaxis = 0.5,xticks = x.ticks,boxsize = 0.1,txt_gp = fpTxtGp(
      ticks = gpar(fontfamily = "", cex = fontsize * 0.8),
      label = gpar(fontfamily = "", cex = fontsize * 0.9),
      summary = gpar(cex = fontsize)
    ),
    colgap = unit(2, "mm"),align = c("l", "l", "l"),mar = unit(c(4,1,9,1), "mm"), #bltr
    clip = c(min(x.ticks), max(x.ticks)), ...
  )
}


#########################################################################################
# CoR Object 1

# forest plots

var.ind=2 # markers location

for (y in c("y1","y2","y3","y4")) {
for (.trt in c(1,0)) {
for (ind in c("b","a")) { # ind a: EIA/PCA use phase 1 data; ind b: EIA/PCA use phase 2 data
#y="y1"; .trt=1; ind="a"; x="EIA.log10d14overd0"
    
    markers=c(t(outer(assays,"."%.%times[c(1,2,4,3)], paste0))); markers; names(markers)=markers    
    fits=list()
    for (x in markers) {
#        if (!contain(x,"cord")) {
#            dat.wide.pop=subset(dat.wide, ppimmfl=="Y")
#            dat.wide.v.pop=subset(dat.wide.v, ppimmfl=="Y")
#            dat.wide.v.pop$wt=dat.wide.v$wt.ppimmfl
#        } else {
#            dat.wide.pop=subset(dat.wide, ppimifl=="Y")
#            dat.wide.v.pop=subset(dat.wide.v, ppimifl=="Y")
#            dat.wide.v.pop$wt=dat.wide.v$wt.ppimifl
#        }
        dat.wide.pop=dat.wide
        dat.wide.v.pop=dat.wide.v
    
        formula=paste0(y, "~", ifelse(contain(x,"log10"),x,"log10("%.%x%.%")"))
        # add maternal risk scores
        if (y=="y1") formula=paste0(formula,"+riskScore.mat.endpoint1") else formula=paste0(formula,"+riskScore.mat.endpoint2")
        # adj vacc2birth for maternal timepoints
        if(!contain(x,"cord")) formula=paste0(formula,"+vacc2birth")
        formula=as.formula(formula); print(formula)
        
        if(ind=="b" | startsWith(x,"RSV")) {
            if (y=="y1") {
                fit=tps.rsv(formula=formula, trt.grp=.trt) # tps.rsv does not work well with many empty cells
            } else {
                dstrat<-svydesign(id=~1,strata=~stratum, weights=~wt, data=subset(dat.wide.v.pop, trt==.trt))
                fit=svyglm(formula, design=dstrat, family="binomial")
            }                                
            #fit=glm(formula, subset(dat.wide.v, trt==.trt), family="binomial", weights=subset(dat.wide.v, trt==.trt)$wt)
        } else {
                fit=glm(formula=formula, subset(dat.wide.pop, trt==.trt), family="binomial")
        }    
        print(summary(fit))
        fits[[x]]=fit
    }
    
    res=c(getFormattedSummary(fits, robust=F, exp=T, rows=var.ind))    
    tab=matrix(res, ncol=length(assays), dimnames=list(paste0(ifelse(contain(times[c(1,2,4,3)],"log10"),"","log10 "),times[c(1,2,4,3)]), assays)); print(tab)
    mytex(tab, file=paste0("tables/CoR_obj1", ind, "_", y, "_trt",.trt), align="c", input.folder=save.results.to)
    
    # multitesting adjustment
    pvals=sapply(fits, function(fit) last(getFixedEf(fit,robust=F)[2,])) # last number in the second row is p value
    tab.2=cbind(pvals, Holm=p.adjust(pvals, method="holm"), BH=p.adjust(pvals, method="BH")); tab.2
    mytex(tab.2, file=paste0("tables/CoR_obj1", ind, "_", y, "_trt",.trt, "_pvals"), align="c", input.folder=save.results.to)
    
    # for some reason the forest plot pdfs are not generated properly when run in a loop, but runs okay when run one by one
    
    # forest plot
    if (ind=="a") {
        mypdf(oma=c(1,0,0,0), onefile=F, width=12,height=6, file=paste0(save.results.to, "/hr_forest", "_", y, "_trt",.trt))  # onefile has to be F otherwise there will be an empty page inserted
            est.ci = sapply(fits[c(which(endsWith(names(fits), "d14")),  which(endsWith(names(fits), "overd0")))] , function (fit) {
                if (length(fit)==1) return (rep(NA,4))
                tmp=getFixedEf(fit, exp=T, robust=F)
                #print(tmp)
                tmp[var.ind,c("OR", "(lower", "upper)", "p.value")]
            })
            colnames(est.ci)=sub(".log10d14overd0"," fold-rise",colnames(est.ci))
            colnames(est.ci)=sub(".d14"," D14",colnames(est.ci))
            theforestplot (point.estimates=est.ci[1,], lower.bounds=est.ci[2,], upper.bounds=est.ci[3,], p.values=est.ci[4,], graphwidth=unit(120, "mm"), table.labels = c("Biomarker", "OR (95% CI)","P-value"), 
                group=colnames(est.ci), decimal.places=2, title=paste0("2-Phase Logistic Regression Adjusting for Baseline Factors: ",ifelse(.trt==1,"Vaccine","Placebo"), " Arm Endpoint ", substr(y,2,100)))
        dev.off()
        
        
        # cord
        mypdf(oma=c(1,0,0,0), onefile=F, width=12,height=6, file=paste0(save.results.to, "/hr_forest_cord", "_", y, "_trt",.trt))  # onefile has to be F otherwise there will be an empty page inserted
            est.ci = sapply(fits[c(which(endsWith(names(fits), "cord")))] , function (fit) {
                if (length(fit)==1) return (rep(NA,4))
                tmp=getFixedEf(fit, exp=T, robust=F)
                #print(tmp)
                tmp[var.ind,c("OR", "(lower", "upper)", "p.value")]
            })
            colnames(est.ci)=sub(".cord"," Cord",colnames(est.ci))
            theforestplot (point.estimates=est.ci[1,], lower.bounds=est.ci[2,], upper.bounds=est.ci[3,], p.values=est.ci[4,], graphwidth=unit(120, "mm"), table.labels = c("Biomarker", "OR (95% CI)","P-value"), 
                group=colnames(est.ci), decimal.places=2, title=paste0("2-Phase Logistic Regression Adjusting for Baseline Factors: ",ifelse(.trt==1,"Vaccine","Placebo"), " Arm Endpoint ", substr(y,2,100)))
        dev.off()
        
    }
    
}
}
}


########################################################################################################
# implement Keith suggestions looking at the diff between cord and baseline 

var.ind=2 # markers location
for (y in c("y2")) {
for (.trt in c(1,0)) {
#y="y2"; .trt=1
ind="a"     
    fits=list()
    for (x in assays) {
#        if (!contain(x,"cord")) {
#            dat.wide.pop=subset(dat.wide, ppimmfl=="Y")
#            dat.wide.v.pop=subset(dat.wide.v, ppimmfl=="Y")
#            dat.wide.v.pop$wt=dat.wide.v$wt.ppimmfl
#        } else {
#            dat.wide.pop=subset(dat.wide, ppimifl=="Y")
#            dat.wide.v.pop=subset(dat.wide.v, ppimifl=="Y")
#            dat.wide.v.pop$wt=dat.wide.v$wt.ppimifl
#        }
        dat.wide.pop=dat.wide
        dat.wide.v.pop=dat.wide.v
        
        dat.wide.pop[["a"]]  =dat.wide.pop[[x%.%".log10cord"]]   - dat.wide.pop[[x%.%".log10d0"]]
        dat.wide.v.pop[["a"]]=dat.wide.v.pop[[x%.%".log10cord"]] - dat.wide.v.pop[[x%.%".log10d0"]]
            
        formula=paste0(y, "~a")        
        if (y=="y1") formula=paste0(formula,"+riskScore.mat.endpoint1") else formula=paste0(formula,"+riskScore.mat.endpoint2") # add maternal risk scores
        # not adj vacc2birth for cord timepoints
        formula=as.formula(formula); print(formula)
        
        if(ind=="b" | startsWith(x,"RSV")) {
            if (y=="y1") {
                #fit=tps.rsv(formula=formula, trt.grp=.trt) # tps.rsv does not work well with many empty cells
            } else {
                dstrat<-svydesign(id=~1,strata=~stratum, weights=~wt, data=subset(dat.wide.v.pop, trt==.trt))
                fit=svyglm(formula, design=dstrat, family="binomial")
            }                                
            #fit=glm(formula, subset(dat.wide.v, trt==.trt), family="binomial", weights=subset(dat.wide.v, trt==.trt)$wt)
        } else {
                fit=glm(formula=formula, subset(dat.wide.pop, trt==.trt), family="binomial")
        }    
        print(summary(fit))
        fits[[x]]=fit
    }
    
    res=c(getFormattedSummary(fits, robust=F, exp=T, rows=var.ind))    
    tab=matrix(res, ncol=length(assays), dimnames=list("log10cordoverd0", assays)); print(tab)
    mytex(tab, file=paste0("tables/CoR_obj1", ind, "_", y, "_trt",.trt, "_cordoverd0"), align="c", input.folder=save.results.to)
    
}
}


########################################################################################################
# scale markers

var.ind=2 # markers location in regression models
.trt=1
y="y2"
tab=NULL
for (t in c("d14","cord")) {
#t="cord"
    
    fits=list()
    for (x in assays) {
        dat.wide.pop=dat.wide
        dat.wide.v.pop=dat.wide.v
        
        dat.wide.pop[["a"]]  =dat.wide.pop[[x%.%".log10"%.%t]]   - dat.wide.pop[[x%.%".log10d0"]]
        dat.wide.v.pop[["a"]]=dat.wide.v.pop[[x%.%".log10"%.%t]] - dat.wide.v.pop[[x%.%".log10d0"]]
            
        formula=paste0(y, "~scale(a)")        
        if (y=="y1") formula=paste0(formula,"+riskScore.mat.endpoint1") else formula=paste0(formula,"+riskScore.mat.endpoint2") # add maternal risk scores
        # adj vacc2birth for maternal timepoints
        if(!contain(t,"cord")) formula=paste0(formula,"+vacc2birth")
        formula=as.formula(formula); print(formula)
        
        if(startsWith(x,"RSV")) {
            if (y=="y1") {
                #fit=tps.rsv(formula=formula, trt.grp=.trt) # tps.rsv does not work well with many empty cells
            } else {
                dstrat<-svydesign(id=~1,strata=~stratum, weights=~wt, data=subset(dat.wide.v.pop, trt==.trt))
                fit=svyglm(formula, design=dstrat, family="binomial")
            }                                
            #fit=glm(formula, subset(dat.wide.v, trt==.trt), family="binomial", weights=subset(dat.wide.v, trt==.trt)$wt)
        } else {
                fit=glm(formula=formula, subset(dat.wide.pop, trt==.trt), family="binomial")
        }    
        print(summary(fit))
        fits[[x]]=fit
    }
    
    res=c(getFormattedSummary(fits, robust=F, exp=T, rows=var.ind))    
    tab=rbind(tab, res)
    
}

colnames(tab)=assays
rownames(tab)=c("log10d14overd0", "log10cordoverd0")
print(tab)
mytex(tab, file=paste0("tables/CoR_obj1_y2_trt1_foldchanges"), align="c", input.folder=save.results.to)






## report all covariates in the regression model, only y2
y="y2"
for (t in c("log10d14overd0","log10d14")) {    
    fits=list()
    for (.trt in c(1,0)) {
        fits[[.trt+1]]=list()
        for (x in assays) {
            formula=paste0(y, "~",x,".",t)
            # add maternal risk scores
            if (y=="y1") formula=paste0(formula,"+riskScore.mat.endpoint1") else formula=paste0(formula,"+riskScore.mat.endpoint2")
            # adj vacc2birth for maternal timepoints
            if(!contain(x,"cord")) formula=paste0(formula,"+vacc2birth")
            formula=as.formula(formula); print(formula)
            
            if(startsWith(x,"RSV")) {
                if (y=="y1") {
                    fit=tps.rsv(formula=formula, trt.grp=.trt)
                } else {
                    dstrat<-svydesign(id=~1,strata=~stratum, weights=~wt.ppimmfl, data=subset(dat.wide.v, trt==.trt & ppimmfl=="Y"))
                    fit=svyglm(formula, design=dstrat, family="binomial")
                }                                
                #fit=glm(formula, subset(dat.wide.v, trt==.trt), family="binomial", weights=subset(dat.wide.v, trt==.trt)$wt)
            } else {
                fit=glm(formula=formula, subset(dat.wide, trt==.trt & ppimmfl=="Y"), family="binomial")
            }    
            fits[[.trt+1]][[x]]=fit
        }
    }
        
    for (.trt in c(1,0)) {
        tab=getFormattedSummary(fits[[1+.trt]], robust=F, exp=T, rows=2:4)        
        rownames(tab)=sub("RSVA.","",rownames(tab))
        mytex(tab, file=paste0("CoR_obj1_", y, "_trt",.trt, "_full_",t), align="c", input.folder=save.results.to)
    }
}



########################################################################################
# marginal risk plots


lwd=1.5
plac.prev=subset(dat.wide, trt==0, y2, drop=T) 
plac.prev=binconf(sum(plac.prev), length(plac.prev))   
ylab="marginal risk of endpoint 2"

mypdf(file=paste0(save.results.to, "/marginal_risk_EIA.log10d14overd0"), width=5, height=5)
    tmp=risks[[1]]
    plot(prob~marker, tmp, type="l", xlab="EIA.log10d14overd0", ylab=ylab, lwd=lwd, ylim=range(c(0,plac.prev, tmp$prob)))
    lines(tmp$marker, tmp$lb, lty=3)
    lines(tmp$marker, tmp$ub, lty=3)
    abline(h=plac.prev, col="gray", lty=c(1,2,2), lwd=lwd)
dev.off()    

mypdf(file=paste0(save.results.to, "/marginal_risk_EIA"), width=5, height=5)
    tmp=cbind(EIA.log10d14overd0=risks[[1]]$prob, EIA.log10d14=risks[[2]]$prob)
    mymatplot(seq(.05,.95,by=0.01), tmp, type="l", legend.x=3, col=1:2, lty=1, xlab="quantile", ylab=ylab, text.width =.4, lwd=lwd, ylim=range(c(0,plac.prev, tmp))) #draw.x.axis=F, xaxt="n"
    abline(h=plac.prev, col="gray", lty=c(1,2,2), lwd=lwd)
    lines(seq(.05,.95,by=0.01), risks[[1]]$lb, lty=3); lines(seq(.05,.95,by=0.01), risks[[1]]$ub, lty=3)
    lines(seq(.05,.95,by=0.01), risks[[2]]$lb, lty=3, col=2); lines(seq(.05,.95,by=0.01), risks[[2]]$ub, lty=3, col=2)
dev.off()    

mypdf(file=paste0(save.results.to, "/marginal_risk_EIA.log10d14overd0_adjd0"), width=5, height=5)
    tmp=cbind("d0 not adjusted"=risks[[1]]$prob, "d0 adjusted"=risks[[3]]$prob)
    mymatplot(risks[[1]]$marker, tmp, type="l", legend.x=3, col=1:2, lty=1, xlab="EIA.log10d14overd0", ylab=ylab, text.width =.5, lwd=lwd, ylim=range(c(0,plac.prev, tmp))) #draw.x.axis=F, xaxt="n"
    abline(h=plac.prev, col="gray", lty=c(1,2,2), lwd=lwd)
    lines(risks[[1]]$marker, risks[[1]]$lb, lty=3); lines(risks[[1]]$marker, risks[[1]]$ub, lty=3)
    lines(risks[[1]]$marker, risks[[3]]$lb, lty=3, col=2); lines(risks[[1]]$marker, risks[[3]]$ub, lty=3, col=2)
dev.off()    
