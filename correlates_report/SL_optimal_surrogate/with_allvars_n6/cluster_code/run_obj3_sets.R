## load required libraries and functions
library("methods")
library("SuperLearner")
library("e1071")
library("glmnet")
library("xgboost")
library("earth")
library("dplyr")
library("conflicted")
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

# code_dir<-'/home/wzhang/Desktop/risk_score/sl_batch_code'
# code_dir <- "/home/bborate/RSVcorr/objective3"
# setwd(code_dir)
num_cores <- parallel::detectCores()
source("sl_screens.R") # set up the screen/algorithm combinations
source("utils.R") # get CV-AUC for all algs
# filter() is used in both dplyr and stats, so need to set the preference to dplyr
conflict_prefer("filter", "dplyr")
conflict_prefer("summarise", "dplyr")


###########################################################################
# impute missing values for covars in maternal and birth variable set
# impute for placebo and vaccine grp separately for risk score analysis
#     use data from all available study subjects in the groups
# dat.covar.imp = read.csv('H:/RSV/dat.covar.imp_v2.0.csv',header = TRUE) 
dat.covar.imp = dat.wide %>%
  mutate(smoker = case_when(smoker == "Y" ~ 1,
                            smoker == "N" ~ 0),
         child5 = case_when(child5 == "Y" ~ 1,
                            child5 == "N" ~ 0),
         daycare = case_when(daycare == "Y" ~ 1,
                             daycare == "N" ~ 0),
         m.ast = case_when(m.ast == "Y" ~ 1,
                           m.ast == "N" ~ 0),
         mhsmopr = case_when(mhsmopr == "Y" ~ 1,
                             mhsmopr == "N" ~ 0))

library(mice)
set.seed(1)
n.imp=1
covars = c("age.at.trt", "bmi", "child5", "smoker", "daycare")
for (.trt in 0:1) {
  imp=mice(dat.covar.imp[dat.covar.imp$trt==.trt, covars],m=n.imp)
  # use the first imputation by default
  dat.covar.imp[dat.covar.imp$trt==.trt, covars]=mice::complete(imp, action=1)
  #
  for (i in 1:n.imp) {
    dat.covar.imp[dat.covar.imp$trt==.trt, covars%.%".imp"%.%i]=mice::complete(imp, action=i)
  }
}

dat.covar.imp <- dat.covar.imp %>%
  mutate(age.at.trt.cat = case_when(age.at.trt > 28 ~ 1,
                                    age.at.trt <= 28 ~ 0))
###########################################################################
# Create combined new dataset which has imputed values of demographics (for phase 1 data) from dat.covar.imp AND 
# imputed values for markers (for phase 2 data) from dat.wide.v
dat.ph1 <- dat.covar.imp %>%
  select(pair.id:stratum, age.at.trt.cat, wt)

dat.ph2 = dat.covar.imp %>%
  filter(sampled=="Y") %>%
  select(pair.id:p.igr, age.at.trt.cat) %>%
  full_join(dat.wide.v %>% select(pair.id, EIA.d0:RSVB.d14, 
                                  EIA.log10d14overd0, PCA.log10d14overd0, RSVA.log10d14overd0, RSVB.log10d14overd0,
                                  sampdy.d0, sampdy.d14, 
                                  EIA.cord:RSVB.cord, 
                                  sampdy.cord:stratum, wt), by = "pair.id") %>%
  mutate(EIA.log10d0 = log10(EIA.d0), EIA.log10d14 = log10(EIA.d14), EIA.log10cord = log10(EIA.cord),
         PCA.log10d0 = log10(PCA.d0), PCA.log10d14 = log10(PCA.d14), PCA.log10cord = log10(PCA.cord),
         RSVA.log10d0 = log10(RSVA.d0), RSVA.log10d14 = log10(RSVA.d14), RSVA.log10cord = log10(RSVA.cord),
         RSVB.log10d0 = log10(RSVB.d0), RSVB.log10d14 = log10(RSVB.d14), RSVB.log10cord = log10(RSVB.cord),
         EIA.2fold.d14overd0 = ifelse(EIA.log10d14 > EIA.log10d0 + log10(2), 1, 0),
         EIA.4fold.d14overd0 = ifelse(EIA.log10d14 > EIA.log10d0 + log10(4), 1, 0),
         PCA.2fold.d14overd0 = ifelse(PCA.log10d14 > PCA.log10d0 + log10(2), 1, 0),
         PCA.4fold.d14overd0 = ifelse(PCA.log10d14 > PCA.log10d0 + log10(4), 1, 0),
         RSVA.2fold.d14overd0 = ifelse(RSVA.log10d14 > RSVA.log10d0 + log10(2), 1, 0),
         RSVA.4fold.d14overd0 = ifelse(RSVA.log10d14 > RSVA.log10d0 + log10(4), 1, 0),
         RSVB.2fold.d14overd0 = ifelse(RSVB.log10d14 > RSVB.log10d0 + log10(2), 1, 0),
         RSVB.4fold.d14overd0 = ifelse(RSVB.log10d14 > RSVB.log10d0 + log10(4), 1, 0))

Z_plus_weights <- dat.ph1 %>% 
  select(pair.id, y1, wt, trt, trt.label, age.at.trt.cat, age.at.trt, bmi, mhsmopr, m.ast, child5, season, smoker, daycare) %>%
  #mutate(ptid = as.numeric(ptid)) %>% 
  filter(!is.na(y1), !is.na(age.at.trt), !is.na(age.at.trt.cat), !is.na(bmi),
         !is.na(mhsmopr), !is.na(m.ast), !is.na(child5), !is.na(season), 
         !is.na(smoker), !is.na(daycare))

# dat.covar.imp.v3 = dat.covar.imp.v2 %>% 
#   rbind(dat.wide.v2)
###########################################################################

# Define treatment: either placebo/vaccine group for objective 3
# Define endpoint: either y1/y2/y3 for objective 3
treatment = 0
endpoint = "y1"
dat.ph1 <- dat.ph1 %>%  filter(trt == treatment)     # 784 records for placebo
dat.ph2 <- dat.ph2 %>%  filter(trt == treatment)     # 155 records for placebo
endpoint.ph1 = (dat.ph1 %>% select(matches(endpoint)))[[1]]
endpoint.ph2 = (dat.ph2 %>% select(matches(endpoint)))[[1]]
if(treatment == 0){
  trt_string = "placebo"
} else {
  trt_string = "vaccine"
}

###########################################################################

matvars <- dat.ph1 %>% 
  select(age.at.trt.cat, age.at.trt, bmi, mhsmopr, m.ast, child5, season, smoker, daycare) %>%
  colnames()

markers <- dat.ph2 %>%
  select(EIA.log10d0, EIA.log10d14, EIA.log10d14overd0, EIA.log10cord, 
         PCA.log10d0, PCA.log10d14, PCA.log10d14overd0, PCA.log10cord, 
         RSVA.log10d0, RSVA.log10d14, RSVA.log10d14overd0, RSVA.log10cord, 
         RSVB.log10d0, RSVB.log10d14, RSVB.log10d14overd0, RSVB.log10cord,
         EIA.2fold.d14overd0, EIA.4fold.d14overd0, 
         PCA.2fold.d14overd0, PCA.4fold.d14overd0, 
         RSVA.2fold.d14overd0, RSVA.4fold.d14overd0, 
         RSVB.2fold.d14overd0, RSVB.4fold.d14overd0) %>%
  # select(-PCA.4fold.d14overd0) %>% # drop PCA.4fold.d14overd0 variable for placebo group as all values are 0!
  colnames()

#####################################################################################################################
## Create variable sets and set up X, Y for super learning
nofold_markers <- markers[1:16]
fold_markers <- markers[17:24]
# Maternal enrollment variables are default in all sets
# 1. None (No markers; only maternal enrollment variables), phase 1 data
varset_none <- rep(FALSE, length(markers))

# 2--17. sets with each individual marker, phase 2 data
for(x in 1:length(nofold_markers)){
  varset <- rep(FALSE, length(markers))
  varset[x] = TRUE
  assign(paste0("varset_", nofold_markers[x]), varset)
}

# 18--21. By time point across assay types, phase 2 data
varset_d0_all4 <- create_varsets(markers, c("EIA.log10d0", "PCA.log10d0", "RSVA.log10d0", "RSVB.log10d0"))
varset_d14_all4 <- create_varsets(markers, c("EIA.log10d14", "PCA.log10d14", "RSVA.log10d14", "RSVB.log10d14"))
varset_d14overd0_all4 <- create_varsets(markers, c("EIA.log10d14overd0", "PCA.log10d14overd0", "RSVA.log10d14overd0", "RSVB.log10d14overd0"))
varset_cord_all4 <- create_varsets(markers, c("EIA.log10cord", "PCA.log10cord", "RSVA.log10cord", "RSVB.log10cord"))

# 22--25. By assay type across timepoints, phase 2 data
varset_EIA_all4 <- create_varsets(markers, c("EIA.log10d0", "EIA.log10d14", "EIA.log10d14overd0", "EIA.log10cord"))
varset_PCA_all4 <- create_varsets(markers, c("PCA.log10d0", "PCA.log10d14", "PCA.log10d14overd0", "PCA.log10cord"))
varset_RSVA_all4 <- create_varsets(markers, c("RSVA.log10d0", "RSVA.log10d14", "RSVA.log10d14overd0", "RSVA.log10cord"))
varset_RSVB_all4 <- create_varsets(markers, c("RSVB.log10d0", "RSVB.log10d14", "RSVB.log10d14overd0", "RSVB.log10cord"))

# 26--28. Maternal baseline markers AND (either d14 or fold-rise or cord markers)
varset_d0_d14 <- create_varsets(markers, c("EIA.log10d0", "PCA.log10d0", "RSVA.log10d0", "RSVB.log10d0", 
                                           "EIA.log10d14", "PCA.log10d14", "RSVA.log10d14", "RSVB.log10d14"))
varset_d0_foldrise <- create_varsets(markers, c("EIA.log10d0", "PCA.log10d0", "RSVA.log10d0", "RSVB.log10d0", 
                                                "EIA.log10d14overd0", "PCA.log10d14overd0", "RSVA.log10d14overd0", "RSVB.log10d14overd0"))
varset_d0_cord <- create_varsets(markers, c("EIA.log10d0", "PCA.log10d0", "RSVA.log10d0", "RSVB.log10d0", 
                                            "EIA.log10cord", "PCA.log10cord", "RSVA.log10cord", "RSVB.log10cord"))

# 29. All markers
varset_all16 <- create_varsets(markers, markers[1:16])

varset_names <- c("1_only_matvars", 
                  "2_varset_EIA.log10d0", "3_varset_EIA.log10d14", "4_varset_EIA.log10d14overd0", "5_varset_EIA.log10cord",
                  "6_varset_PCA.log10d0", "7_varset_PCA.log10d14", "8_varset_PCA.log10d14overd0", "9_varset_PCA.log10cord", 
                  "10_varset_RSVA.log10d0", "11_varset_RSVA.log10d14", "12_varset_RSVA.log10d14overd0", "13_varset_RSVA.log10cord",
                  "14_varset_RSVB.log10d0", "15_varset_RSVB.log10d14", "16_varset_RSVB.log10d14overd0", "17_varset_RSVB.log10cord", 
                  "18_varset_d0_all4","19_varset_d14_all4", "20_varset_d14overd0_all4", "21_varset_cord_all4",
                  "22_varset_EIA_all4", "23_varset_PCA_all4", "24_varset_RSVA_all4", "25_varset_RSVB_all4",
                  "26_varset_d0_d14", "27_varset_d0_foldrise", "28_varset_d0_cord",
                  "29_varset_all16")

## set up a matrix of all 
varset_matrix <- rbind(varset_none, 
                       varset_EIA.log10d0, varset_EIA.log10d14, varset_EIA.log10d14overd0, varset_EIA.log10cord,
                       varset_PCA.log10d0, varset_PCA.log10d14, varset_PCA.log10d14overd0, varset_PCA.log10cord,
                       varset_RSVA.log10d0, varset_RSVA.log10d14, varset_RSVA.log10d14overd0, varset_RSVA.log10cord, 
                       varset_RSVB.log10d0, varset_RSVB.log10d14, varset_RSVB.log10d14overd0, varset_RSVB.log10cord,
                       varset_d0_all4, varset_d14_all4, varset_d14overd0_all4, varset_cord_all4,
                       varset_EIA_all4, varset_PCA_all4, varset_RSVA_all4, varset_RSVB_all4,
                       varset_d0_d14, varset_d0_foldrise, varset_d0_cord,
                       varset_all16)

job_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
this_var_set <- varset_matrix[job_id, ]
cat("\n Running ", varset_names[job_id], "\n")

#############################################################################################################
X_covars2adjust_ph1 <- dat.ph1 %>% select(all_of(matvars))
X_covars2adjust_ph2 <- dat.ph2 %>% select(all_of(c(matvars, markers)))   

## scale covars2adjust to have mean 0, sd 1 for all vars
for (a in colnames(X_covars2adjust_ph1)) {
  X_covars2adjust_ph1[[a]] <- scale(X_covars2adjust_ph1[[a]],
                                    center = mean(X_covars2adjust_ph1[[a]], na.rm = T), 
                                    scale = sd(X_covars2adjust_ph1[[a]], na.rm = T))    
}


for (a in colnames(X_covars2adjust_ph2)) {
  X_covars2adjust_ph2[[a]] <- scale(X_covars2adjust_ph2[[a]],
                                    center = mean(X_covars2adjust_ph2[[a]], na.rm = T), 
                                    scale = sd(X_covars2adjust_ph2[[a]], na.rm = T))    
}


X_covars2adjust_ph2 <- X_covars2adjust_ph2 %>% 
  select_if(function(x) any(!is.na(x))) # if placebo group and endpoint is y2, PCA.4fold.d14overd0 variable
                                        # has 0 variance, and returns all NAN's from scale function. This drops such variables.

##############################################################################################################
# select data based on job_id
# if(job_id == 1){
#   X_markers_varset <- X_covars2adjust_ph2[1:9]
#   Y = endpoint.ph2
#   weights = dat.ph2$wt
#   # X_markers_varset <- X_covars2adjust_ph1
#   # Y = endpoint.ph1
#   # weights = rep(1, length(Y))
#   sl_lib <- c("SL.mean", "SL.glm","SL.glm.interaction","SL.bayesglm", "SL.step", "SL.glmnet","SL.gam","SL.cforest","SL.xgboost")
#   #sl_lib <- SL_library
# } else {
  X_markers_varset <- bind_cols(X_covars2adjust_ph2[1:9], 
                                X_covars2adjust_ph2[10:length(X_covars2adjust_ph2)] %>% select(names(X_covars2adjust_ph2[10:25])[this_var_set]),
                                X_covars2adjust_ph2[26:length(X_covars2adjust_ph2)])
  Y = endpoint.ph2
  weights = dat.ph2$wt
  sl_lib <- SL_library
# }

# C <- rep(1, length(Y))
# Z <- data.frame(Y = Y, X_markers_varset %>% select(age.at.trt, age.at.trt.cat, bmi, mhsmopr, m.ast, child5, season, smoker, daycare))

treatmentDAT <- dat.ph2 %>% select(pair.id, trt, wt, y1, all_of(c(matvars, markers))) %>%
  filter(trt == treatment) %>%
  select(-trt)

# match the rows in treatmentDAT to get Z, C
all_cc_treatment <- Z_plus_weights %>%
  filter(pair.id %in% treatmentDAT$pair.id, trt == treatment)
# pull out the participants who are NOT in the cc cohort and received the vaccine
all_non_cc_treatment <- Z_plus_weights %>%
  filter(!(pair.id %in% treatmentDAT$pair.id), trt == treatment)
# put them back together
phase_1_data_treatmentDAT <- dplyr::bind_rows(all_cc_treatment, all_non_cc_treatment) %>%
  select(-trt, -trt.label)
Z_treatmentDAT <- phase_1_data_treatmentDAT %>%
  select(-pair.id, -wt)
all_ipw_weights_treatment <- phase_1_data_treatmentDAT %>%
  pull(wt)
C <- (phase_1_data_treatmentDAT$pair.id %in% treatmentDAT$pair.id)



# else if (job_id !=1 & job_id < 26) {
#   X_markers_varset <- bind_cols(X_covars2adjust_ph2[1:9], 
#                                 X_covars2adjust_ph2[10:length(X_covars2adjust_ph2)] %>% select(names(X_covars2adjust_ph2[10:25])[this_var_set]),
#                                 X_covars2adjust_ph2[26:length(X_covars2adjust_ph2)])
#   Y = endpoint.ph2
#   weights = dat.ph2$wt
#   sl_lib <- c("SL.mean", "SL.glm","SL.glm.interaction","SL.bayesglm", "SL.step", "SL.glmnet","SL.gam","SL.cforest","SL.xgboost")
# } 
##############################################################################################################

## set up outer folds for cv variable importance; do stratified sampling
V_outer <- 5
V_inner <- length(Y) - 1

## ---------------------------------------------------------------------------------
## run super learner, with leave-one-out cross-validation and all screens
## do 10 random starts, average over these
## use assay groups as screens
## ---------------------------------------------------------------------------------
## ensure reproducibility
set.seed(20201202)
# seeds <- round(runif(1, 1000, 10000)) # do only one seed as trial
seeds <- round(runif(10, 1000, 10000)) # average over 10 random starts


##solve cores issue
library(RhpcBLASctl)
blas_get_num_procs()
blas_set_num_threads(1)
print(blas_get_num_procs())
stopifnot(blas_get_num_procs()==1)

fits <- parallel::mclapply(seeds, FUN = run_cv_sl_once, 
                           Y = Y, 
                           X_mat = X_markers_varset, 
                           family = "binomial",
                           obsWeights = weights,
                           all_weights = all_ipw_weights_treatment,
                           sl_lib = sl_lib,
                           method = "method.CC_nloglik",
                           cvControl = list(V = V_outer, stratifyCV = TRUE),
                           innerCvControl = list(list(V = V_inner)),
                           Z = Z_treatmentDAT, 
                           C = C, 
                           # z_lib = c("SL.glm", "SL.bayesglm", "SL.step", "SL.gam","SL.cforest"), # new arguments
                           z_lib = "SL.glm",
                           scale = "logit", # new argument
                           vimp = FALSE,
                           mc.cores = num_cores
)

saveRDS(fits, paste0("/fh/fast/gilbert_p/RSV/objective3/results/phase2_with_allvars_n6/", endpoint, "_", trt_string, "/slfits_", endpoint, "_", trt_string, "_", varset_names[job_id], ".rds"))




# X_mat = X_markers_varset
# family = "binomial"
# Z = Z_treatmentDAT
# z_lib = c("SL.glm", "SL.bayesglm", "SL.step", "SL.gam","SL.cforest")
# z_lib = "SL.glm"
# obsWeights = weights
# all_weights = all_ipw_weights_treatment
# scale = "logit"
# method = "method.CC_nloglik"
# cvControl = list(V = V_outer, stratifyCV = TRUE)
# innerCvControl = list(list(V = V_inner))
# vimp = FALSE
# mc.cores = num_cores