# ---------------------------------------------------------------------------------
# pre-process the data
# ---------------------------------------------------------------------------------
# read in the full dataset
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

# Consider only placebo data for risk score analysis
# placebo group
dat.wide.plbo <- dat.covar.imp %>% filter(trt == 0)
X_mat_plbo <- dat.wide.plbo  %>%
  select(age.at.trt.cat, age.at.trt, 
         bmi, mhsmopr, m.ast,child5,season, smoker,daycare)  # prebrth (high number of missing values) and hiv (all non-NA values are negative) are not considered !!!!
X_ped_plbo <- dat.wide.plbo  %>%select( p.sex, iwt, p.lbw, 
                   iwtlen, hdcirc, ga, p.small, p.igr)

Y1_plbo <- dat.wide.plbo$y1
Y2_plbo <- dat.wide.plbo$y2 

mat_plbo_id=dat.wide.plbo$pair.id

#vaccine group
dat.wide.vacc<-dat.covar.imp %>% filter(trt == 1)
X_mat_vacc <- dat.wide.vacc  %>%
  select(age.at.trt.cat, age.at.trt, 
         bmi, mhsmopr, m.ast,child5,season, smoker,daycare)  # prebrth (high number of missing values) and hiv (all non-NA values are negative) are not considered !!!!
X_ped_vacc <- dat.wide.vacc  %>%select( p.sex, iwt, p.lbw, 
                 iwtlen, hdcirc, ga, p.small, p.igr)
Y1_vacc <- dat.wide.vacc$y1
Y2_vacc <- dat.wide.vacc$y2  

mat_vacc_id=dat.wide.vacc$pair.id

weights_placebo = rep(1, length(Y1_plbo))
weights_vaccine = rep(1, length(Y1_vacc))
## set up outer folds for cv variable importance; do stratified sampling
V_outer <- 5
V_inner <- length(Y1_plbo) - 1


## construct risk score-----------------------

risk_score_mat_y1=SuperLearner::SuperLearner(Y = Y1_plbo, X = X_mat_plbo, family = "binomial",
                           obsWeights = weights_placebo, SL.library = SL_library_1,
                           method = "method.CC_nloglik", cvControl = list(list(V = V_inner))
                              )


risk_score_mat_y2=SuperLearner::SuperLearner(Y = Y2_plbo, X = X_mat_plbo, family = "binomial",
                           obsWeights = weights_placebo, SL.library = SL_library_2,
                           method = "method.CC_nloglik", cvControl = list(list(V = V_inner)))

risk_score_ped_y1=SuperLearner::SuperLearner(Y = Y1, X = X_ped, family = "binomial",
                           obsWeights = weights, SL.library = SL_library_3,
                           method = "method.CC_nloglik", cvControl = list(list(V = V_inner)))

risk_score_ped_y2=SuperLearner::SuperLearner(Y = Y2, X = X_ped, family = "binomial",
                           obsWeights = weights, SL.library =SL_library_1,
                           method = "method.CC_nloglik", cvControl = list(list(V = V_inner)))

#scores for vaccine in mat and ped
mat_vaccine_y1=predict(risk_score_mat_y1, newdata=X_mat_vacc , onlySL = TRUE)[[1]]
mat_vaccine_y2=predict(risk_score_mat_y2, newdata=X_mat_vacc , onlySL = TRUE)[[1]]

score_ped_y1=predict(risk_score_ped_y1, newdata=X_ped, onlySL = TRUE)
score_ped_y2=predict(risk_score_ped_y2, newdata=X_ped, onlySL = TRUE)



#setwd('/home/wzhang/Desktop/risk_score')
#save RDS
saveRDS(score_mat_y1[[1]],'risk_score_mat_y1.rds')
saveRDS(score_mat_y2[[1]],'risk_score_mat_y2.rds')
saveRDS(score_ped_y1[[1]],'risk_score_ped_y1.rds')
saveRDS(score_ped_y2[[1]],'risk_score_ped_y2.rds')

                                  





## save mat_scores in a dataframe with id---------------------------

## read rds for mat y1 and y2
mat_y1_1=readRDS('/fh/fast/gilbert_p/RSV/risk_score/superlearner/results/maternal_add3/endpoint1/sl_fits_mat1_y1_add3.rds')
mat_y2_2=readRDS('/fh/fast/gilbert_p/RSV/risk_score/superlearner/results/maternal_add3/endpoint2/sl_fits_mat2_y2_add3.rds')

#find the best auc 
mat_max_y1=which.max(sapply(1:10,function(x) mat_y1_1[[x]]$aucs[1,3]))
mat_max_y2=which.max(sapply(1:10,function(x) mat_y2_2[[x]]$aucs[1,3]))

# scores for placebo in mat
mat_placebo_y1=mat_y1_1[[mat_max_y1]]$fit$SL.predict
mat_placebo_y2=mat_y2_2[[mat_max_y2]]$fit$SL.predict

nvacc=length(mat_vacc_id)
nplbo=length(mat_plbo_id)



mat_allscores=as.data.frame(rbind(cbind(id=mat_vacc_id,score=mat_vaccine_y1,group=rep('vaccine',nvacc),endpoint=rep('y1',nvacc)),
                       cbind(id=mat_vacc_id,score=mat_vaccine_y2,group=rep('vaccine',nvacc),endpoint=rep('y2',nvacc)),
                       cbind(id=mat_plbo_id,score=mat_placebo_y1,group=rep('placebo',nplbo),endpoint=rep('y1',nplbo)),
                       cbind(id=mat_plbo_id,score=mat_placebo_y2,group=rep('placebo',nplbo),endpoint=rep('y2',nplbo))
                       ))
saveRDS(mat_allscores,'/fh/fast/gilbert_p/RSV/risk_score/mat_scores.rds')
