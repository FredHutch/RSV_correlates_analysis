# Purpose: Principal surrogate evaluation (CoP objective 2 in the RSV correlates SAP)
#          This code computes point estimates of the VE curve.
# Method:  Juraska, Huang, and Gilbert (Biostatistics, 2020)
#          R package pssmooth, version 1.0.2
# Input:   R data package RSVcorr
# Output:  RData files with point estimates of the VE curve
# Author:  Michal Juraska

rm(list=ls(all=TRUE))

codeDir <- "h:/SCHARP/RSV/RSVcorrelatesAnalysis/correlates_report"
outDir <- "h:/SCHARP/RSV/RSVcorrelatesAnalysis/correlates_report/Routput"

library(RSVcorr)
# library(pssmooth)
# 'pssmoothForRSV' is a modification of the package code tailored for this analysis
source(file.path(codeDir, "pssmoothForRSV.R"))

# repeat the analysis for the two PP-IMM cohorts
for (cohort in c("PP-IMM-M", "PP-IMM-I")){
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
  
  # endpoints 1 and 2
  for (y in c("y1", "y2")){
    # EIA d14 and fold-change, and PCA d14 and fold-change (all on log10 scale)
    for (ps in c("EIA.log10d14", "EIA.log10d14overd0", "PCA.log10d14", "PCA.log10d14overd0")){
      
      formula <- as.formula(paste0(y, " ~ ", ps, " + factor(vacc2birthLESSTHAN30) + factor(riskScoreBin.mat.endpoint", substr(y, 2, 2), ")"))
      print(formula)
      
      riskCurve(formula, bsm=paste0(substr(ps, 1, 9), "d0"), tx="trt", data=dat.wide.ppim, weights="wt", 
                saveFile=paste0("VEcurve_", cohort, "_", y, "_", ps, ".RData"), saveDir=outDir)
    }
  }
}
