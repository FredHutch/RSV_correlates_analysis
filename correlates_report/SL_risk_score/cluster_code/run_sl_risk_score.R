## load required libraries and functions
library("methods")
library("SuperLearner")
library("e1071")
library("glmnet")
library("xgboost")
library("earth")
library("dplyr")
## only run this if something has changed
# devtools::install_local("RSVcorr.tar.gz")
library("RSVcorr")
library("kyotil")
library("argparse")
library(vimp)
library(nloptr)
library(RhpcBLASctl)

## set up code directory
# if (!is.na(Sys.getenv("RSTUDIO", unset = NA))) { # if running locally
#   code_dir <- "code/"
# } else {
#   code_dir <- ""
# }

code_dir<-'/home/wzhang/Desktop/risk_score/sl_batch_code'
setwd(code_dir)
num_cores <- parallel::detectCores()
source("sl_screens.R") # set up the screen/algorithm combinations
source("utils.R") # get CV-AUC for all algs




# ---------------------------------------------------------------------------------
# pre-process the data
# ---------------------------------------------------------------------------------
# read in the full dataset
#data("dat.wide", package = "RSVcorr")
## read in the super learner variables
## Check if any covariates have corrleation > 0.9
# 
# dat.wide <- dat.wide %>% 
#   mutate(age.at.trt.cat = case_when(age.at.trt > 28 ~ 1,
#                                     age.at.trt <= 28 ~ 0),          # Add this categorical variable to dat.wide and dat.long datasets !!!!
#          mhsmopr = case_when(mhsmopr == "Y" ~ 1,
#                              mhsmopr == "N" ~ 0),
#          m.ast = case_when(m.ast == "Y" ~ 1,
#                            m.ast == "N" ~ 0)) %>%
#   mutate(smoker = case_when(smoker == "Y" ~ 1,
#                             smoker == "N" ~ 0),
#          p.sex = case_when(p.sex == "F" ~ 1,
#                            p.sex == "M" ~ 0),
#          child5 = case_when(child5 == "Y" ~ 1,
#                             child5 == "N" ~ 0),
#          daycare = case_when(daycare == "Y" ~ 1,
#                              daycare == "N" ~ 0),
#          p.lbw = case_when(p.lbw == "Y" ~ 1,
#                            p.lbw == "N" ~ 0),
#          ga = ga * 7,                                      # Convert ga to days in data package !!!!
#          p.small = case_when(p.small == "Y" ~ 1,
#                              p.small == "N" ~ 0),
#          p.igr = case_when(p.igr == "Y" ~ 1,
#                            p.igr == "N" ~ 0),
#          vacc2birthMORETHAN30 = case_when(vacc2birth >= 30 ~ 1,
#                                           vacc2birth < 30 ~ 0))
# 

###########################################################################
# impute missing values for covars in maternal and birth variable set
# impute for placebo and vaccine grp separately for risk score analysis
#     use data from all available study subjects in the groups
dat.covar.imp=read.csv('dat.covar.imp_v2.0.csv',header = TRUE)
#dat.covar.imp = dat.wide
library(mice)
set.seed(1)
n.imp=1
covars = c("bmi", "smoker", "child5", "daycare", "iwt", "iwtlen", "hdcirc", "ga") 
for (.trt in 0:1) {
  imp=mice(dat.covar.imp[dat.covar.imp$trt==.trt, covars],m=n.imp)
  # use the first imputation by default
  dat.covar.imp[dat.covar.imp$trt==.trt, covars]=mice::complete(imp, action=1)
  #
  for (i in 1:n.imp) {
    dat.covar.imp[dat.covar.imp$trt==.trt, covars%.%".imp"%.%i]=mice::complete(imp, action=i)
  }
}
###########################################################################
# Consider only placebo data for risk score analysis
dat.wide.plbo <- dat.covar.imp %>% filter(trt == 0)
X_mom <- dat.wide.plbo %>%
  select(age.at.trt.cat, age.at.trt, bmi, mhsmopr, m.ast)   # prebrth (high number of missing values) and hiv (all non-NA values are negative) are not considered !!!!
X_ped <- dat.wide.plbo %>%
  select(season, smoker, p.sex, child5, daycare, iwt, p.lbw, iwtlen, hdcirc,
         ga, p.small, p.igr) 
# Should vacc2birthLESSTHAN30 be coded 1 for >=30 days, and 0 for 14-29 days?
X_covars2adjust <- data.frame(X_mom, X_ped)
## scale covars2adjust to have mean 0, sd 1 for all vars
for (a in colnames(X_covars2adjust)) {
  X_covars2adjust[[a]] <- scale(X_covars2adjust[[a]],
                                center = mean(X_covars2adjust[[a]], na.rm = T),  
                                scale = sd(X_covars2adjust[[a]], na.rm = T))    
}
# Taking the composite endpoint as outcome
Y1 <- dat.wide.plbo$y1
Y2 <- dat.wide.plbo$y2                    
placebos <- cbind.data.frame(Y1,Y2, X_covars2adjust) 
# prep for running SL
Y_placebo <- placebos$Y1
X_placebo_mat <- placebos %>% 
  select(age.at.trt.cat, age.at.trt, 
         bmi, mhsmopr, m.ast,child5,season, smoker,daycare)
X_placebo_ped <- placebos %>% 
  select( p.sex, iwt, p.lbw, 
         iwtlen, hdcirc, ga, p.small, p.igr)
weights_placebo = rep(1, length(Y_placebo))
## set up outer folds for cv variable importance; do stratified sampling
V_outer <- 5
V_inner <- length(Y_placebo) - 1


## ---------------------------------------------------------------------------------
## run super learner, with leave-one-out cross-validation and all screens
## do 10 random starts, average over these
## use assay groups as screens
## ---------------------------------------------------------------------------------
## ensure reproducibility
set.seed(4747)
seeds <- round(runif(10, 1000, 10000)) # average over 10 random starts
# mat 1:8, ped 1:8
# # running with high corr scan screen
# fits <- parallel::mclapply(seeds, FUN = run_cv_sl_once, Y = Y_placebo, X = X_placebo, family = "binomial",
#               obsWeights = weights_placebo,
#               sl_lib = SL_library[8],
#               method = "method.CC_nloglik",
#               cvControl = list(V = V_outer, stratifyCV = TRUE),
#               innerCvControl = list(list(V = V_inner)),
#               vimp = FALSE,
#               mc.cores = num_cores
# )
# sl_fits_varset_11_all_1 <- fits

# running with all screens
job_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
cat("\n Running ", job_id, "\n")




##solve cores issue
library(RhpcBLASctl)
blas_get_num_procs()
blas_set_num_threads(1)
print(blas_get_num_procs())
stopifnot(blas_get_num_procs()==1)


#---------------------------------------
# description of job assignment 

#1:8 mat y1
#8:16 mat y2
#17:24 ped y1
#24:32 ped y2



if(job_id<=8){
fits1 <- parallel::mclapply(seeds, FUN = run_cv_sl_once, Y = Y1, X_mat = X_placebo_mat, family = "binomial",
                                                                obsWeights = weights_placebo,
                                                                sl_lib = eval(parse(text=paste0('SL_library_',job_id))),
                                                                method = "method.CC_nloglik",
                                                                cvControl = list(V = V_outer, stratifyCV = TRUE),
                                                                innerCvControl = list(list(V = V_inner)),
                                                                vimp = FALSE,
                                                                mc.cores = num_cores
)
saveRDS(fits1, paste0("/fh/fast/gilbert_p/RSV/risk_score/superlearner/results/maternal_add3/endpoint1/sl_fits_mat", job_id, "_y1_add3.rds"))
} else if(job_id>8 && job_id<=16){
  k=job_id-8
  fits1 <- parallel::mclapply(seeds, FUN = run_cv_sl_once, Y = Y2, X_mat = X_placebo_mat, family = "binomial",
                              obsWeights = weights_placebo,
                              sl_lib = eval(parse(text=paste0('SL_library_',k))),
                              method = "method.CC_nloglik",
                              cvControl = list(V = V_outer, stratifyCV = TRUE),
                              innerCvControl = list(list(V = V_inner)),
                              vimp = FALSE,
                              mc.cores = num_cores
  )
  saveRDS(fits1, paste0("/fh/fast/gilbert_p/RSV/risk_score/superlearner/results/maternal_add3/endpoint2/sl_fits_mat", k, "_y2_add3.rds"))
} else if(job_id>16 && job_id<=24){
  k=job_id-16
  fits1 <- parallel::mclapply(seeds, FUN = run_cv_sl_once, Y = Y1, X_mat = X_placebo_ped, family = "binomial",
                              obsWeights = weights_placebo,
                              sl_lib = eval(parse(text=paste0('SL_library_',k))),
                              method = "method.CC_nloglik",
                              cvControl = list(V = V_outer, stratifyCV = TRUE),
                              innerCvControl = list(list(V = V_inner)),
                              vimp = FALSE,
                              mc.cores = num_cores
  )
  saveRDS(fits1, paste0("/fh/fast/gilbert_p/RSV/risk_score/superlearner/results/pediatric_add3/endpoint1/sl_fits_ped", k, "_y1_add3.rds"))
} else{
  k=job_id-24
  fits1 <- parallel::mclapply(seeds, FUN = run_cv_sl_once, Y = Y2, X_mat = X_placebo_ped, family = "binomial",
                              obsWeights = weights_placebo,
                              sl_lib = eval(parse(text=paste0('SL_library_',k))),
                              method = "method.CC_nloglik",
                              cvControl = list(V = V_outer, stratifyCV = TRUE),
                              innerCvControl = list(list(V = V_inner)),
                              vimp = FALSE,
                              mc.cores = num_cores
  )
  saveRDS(fits1, paste0("/fh/fast/gilbert_p/RSV/risk_score/superlearner/results/pediatric_add3/endpoint2/sl_fits_ped", k, "_y2_add3.rds"))
}


