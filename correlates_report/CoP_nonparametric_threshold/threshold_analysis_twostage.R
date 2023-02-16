library(dplyr)
setwd("~/RSVcorrelatesAnalysis/correlates_report/nonparam_threshold")

#renv::install("~/RSVcorr/RSVcorr_2020.10.31.tar.gz")
renv::install("mvtnorm")
library(RSVcorr); stopifnot(packageVersion("RSVcorr")>="2020.10.31")
source(file.path("tmleThresh.R"))

library(SuperLearner)
#devtools::install_github("tlverse/sl3")
#renv::install("RSVcorr_2020.10.31.tar.gz")
library(sl3)
library(here)
library(data.table)
library(future)
library(kableExtra)

 ##### Change the below for different 
endpoint = "y2" #y1 y2 or y3
treatment = "Vaccine" # "Pooled", "Vaccine", "Placebo"

#### This should take about 10-20 minutes to run. It can be made much faster by  using a single machine learning algorithm as opposed to SuperLearner
#### Options:

covariate_adjusted <- T #### Estimate threshold-response function with covariate adjustment
fast_analysis <- TRUE  ### Perform a fast analysis using glmnet
include_interactions <- FALSE #### Include algorithms that model interactions between covariates
threshold_grid_size <- 25   ### Number of thresholds to estimate (equally spaced in quantiles). Should be 15 at least for the plots of the threshold-response and its inverse to be representative of the true functions.
risks_to_estimate_thresh_of_protection <- NULL  ### Risk values at which to estimate threshold of protection...
###                Leaving at NULL (default) is recommended to ensure risks fall in range of estimate. Example: c(0, 0.0005, 0.001, 0.002, 0.003)

# At least 15 thresholds
threshold_grid_size <- max(threshold_grid_size,15)
###### Running threshold analysis
######

#?RSVcorr::dat.wide
data("dat.wide", package = "RSVcorr")
data("dat.wide.v", package = "RSVcorr")
dat.wide_red <- dat.wide[, c("pair.id", "sampled", "ppintersectfl", "trt"), drop = F]
data <- merge(dat.wide_red, dat.wide.v, all.x = T)



# 12 biomarkers
# here EIA and PCA should have no weights; RSVA and RSVB belong to phase 2, so should have weights
EIA<-colnames(data)[grepl('EIA.log10',colnames(data))][c(2,3,4)]
PCA<-colnames(data)[grepl('PCA.log10',colnames(data))][c(2,3,4)]
RSVA<-colnames(data)[grepl('RSVA.log10',colnames(data))][c(2,3,4)]
RSVB<-colnames(data)[grepl('RSVB.log10',colnames(data))][c(2,3,4)]

immue_variable<-c(EIA,PCA,RSVA,RSVB)

# chosen covariates
covariates<-c("m.ageCAT","age.at.trt","bmi","mhsmopr","m.ast","child5","season","smoker","daycare","vacc2birth")
#data<-data[complete.cases(data[,c(covariates)]),]

data$m.ageCAT<-ifelse(data$m.ageCAT=="> 28 yrs",1,0)
data$mhsmopr<-ifelse(data$mhsmopr=="Y",1,0)
data$m.ast<-ifelse(data$m.ast=="Y",1,0)
data$child5<-ifelse(data$child5=="Y",1,0)
data$smoker<-ifelse(data$smoker=="Y",1,0)
data$daycare<-ifelse(data$daycare=="Y",1,0)

# vaccine group
vaccine<-data[data$trt==1,]

# poled dataset and vacc dataset for EIA and PCA
all<-data[,c(covariates,immue_variable,'y1','y2')]
vacc<-vaccine[,c(covariates,immue_variable,'y1','y2')]


# poled dataset and vacc dataset for RSVA and RSVB
all_wt<-data[,c(covariates,immue_variable,'y1','y2')]
vacc_wt<-vaccine[,c(covariates,immue_variable,'y1','y2')]



# Then use defined "covariates" and "immue_variable" to build model

#--------------------------------------------------------------------------
# For log10 fold-change from maternal baseline to infant cord blood
data<-mutate(data,EIA.log10cordoverd0=EIA.log10cord-EIA.log10d0)
data<-mutate(data,PCA.log10cordoverd0=PCA.log10cord-PCA.log10d0)
data<-mutate(data,RSVA.log10cordoverd0=RSVA.log10cord-RSVA.log10d0)
data<-mutate(data,RSVB.log10cordoverd0=RSVB.log10cord-RSVB.log10d0)
data <- data[ data$ppintersectfl=="Y",]

####################################################
########## Run threshold analysis for raw results
####################################################





library(SuperLearner)
#devtools::install_github("tlverse/sl3")
library(sl3)
library(here)
library(data.table)
library(future)
library(kableExtra)


source(file.path("tmleThresh.R"))

####################################################
#### Data set up
####################################################
# Store the two stage sampling grouping variable at grp
# Reference time for analysis (170 days default)

# Extract outcome (Y) and censoring (Delta) variable

save_file_key <- paste0(endpoint,"_", treatment)
data$outcome <- data[[endpoint]]
folder <- paste0("input/", save_file_key)
if(treatment == "Pooled") {
  data <- data
} else if(treatment == "Vaccine") {
  data <- data[data$trt==1,]
} else if(treatment == "Placebo") {
  data <- data[data$trt==0,]
} else {
  stop("Not valid")
}

data$sampled <- as.numeric(data$sampled=="Y")
# Biomarkers to threshold search'
assays = c("RSVB.log10cordoverd0","RSVA.log10cordoverd0")
assays = c("RSVA.log10cordoverd0","RSVB.log10cordoverd0", "RSVA.log10d14overd0", "RSVB.log10d14overd0", "RSVA.log10cord", "RSVB.log10cord", "RSVA.log10d0", "RSVB.log10d0", "RSVA.log10d14", "RSVB.log10d14") # "RSVB.log10cordoverd0""RSVA.log10cordoverd0")
#assays = c("RSVA.log10d0", "RSVB.log10d0") # "RSVB.log10cordoverd0""RSVA.log10cordoverd0")

#assays <- assays[1]

# Subset to population of interest

# Baseline variables to adjust for
#baseline <- c("Hispanic", "MinorityInd", "HighRiskInd", "BRiskScore", "age.geq.65", "Sex", "race", "BMI") # Add "age"?
baseline <-c("m.ageCAT","age.at.trt","bmi","mhsmopr","m.ast","child5","season","smoker","daycare","vacc2birth")

# Subset to population of interest
#keep <-   data$Trt ==1 & data$Bserostatus==0 & data$Perprotocol==1

# Baseline variables to adjust for
#baseline <- c("Hispanic", "MinorityInd", "HighRiskInd", "BRiskScore", "age.geq.65", "Sex", "race", "BMI") # Add "age"?

####################################################
########## Run threshold analysis for raw results
####################################################
output_list <- list()
# For each marker in assays, run the threshold analysis and save results as csv
for(marker in assays) {

  # Convert marker to marker name found in data
  #marker <- paste0("Day57", marker)
  outcome <- "outcome" #"EventInd"
  #Delta = "Delta"

  print(marker)
  # Subset to necessary variables
  data_full <- data[, c(outcome,  marker, baseline, "wt.ppintersectfl", "sampled", "trt")]
  data_full$wt <- data_full$wt.ppintersectfl
  #keep <- !(is.na(data_full[[marker]]) & data_full$sampled==1)
  #data_full <- data_full[keep,]
  # Subset to those selected for two stage sampling
  data_biased <- data_full[data_full$sampled==1,]

  ####################################################
  #### Machine learning algoriths to use for estimation
  ####################################################
  #glm with bayesian penalization
  bayesglm_sl_lrnr <- make_learner(Lrnr_pkg_SuperLearner, "SL.bayesglm", outcome_type = variable_type("binomial"))
  # glm inluding two way interactions
  lrnr_SL.inter<-  make_learner(Lrnr_pkg_SuperLearner, "SL.glm.interaction", outcome_type = variable_type("binomial"))
  # main term glm
  lrnr_glm <-  Lrnr_glm$new()
  # main term lasso glm
  lrnr_glmnet <- Lrnr_glmnet$new()
  # Empirical mean estimator
  lrnr_mean <- Lrnr_mean$new()
  #baseline <- c("MinorityInd", "HighRiskInd", "Age", "BRiskScore", "age.geq.65")
  interactions <- unlist(apply(combn(baseline,2),2, function(v) list(v)), recursive = F)
  # Lasso with two way interactions
  lrnr_glmnet_inter <- Pipeline$new(Lrnr_define_interactions$new(
    interactions
  ),lrnr_glmnet)
  # bayes glm with interactions
  lrnr_bayes_inter <- Pipeline$new(Lrnr_define_interactions$new(
    interactions
  ),bayesglm_sl_lrnr)

  # SuperLearner of estimators where best is selected with cross validation
  # If too few end points, remove lrnr_bayes_inter and lrnr_glmnet_inter
  # To add glm and and glm with interactions, add lrnr_glm and lrnr_SL.inter

  if(include_interactions) {
    sl_list <- list(lrnr_bayes_inter, lrnr_glmnet, lrnr_mean, bayesglm_sl_lrnr, lrnr_glmnet_inter, lrnr_glm)
  } else {
    sl_list <- list(lrnr_glm, lrnr_glmnet, lrnr_mean, bayesglm_sl_lrnr)
  }
  lrnr <- Lrnr_sl$new(sl_list, Lrnr_cv_selector$new(loss_loglik_binomial))
  # If something goes wrong then use
  # lrnr <- lrnr_glm or
  #lrnr <- lrnr_glmnet
  # For fast results do
  # lrnr <- lrnr_glmnet_inter
  # For unadjusted results
  #lrnr <- Lrnr_mean$new()

  if(fast_analysis) {
    if(include_interactions) {
      lrnr <- lrnr_glmnet_inter
    } else {
      lrnr <- lrnr_glmnet
    }

  }

  if(!covariate_adjusted) {
    lrnr <- Lrnr_mean$new()
  }

  ####################################################
  ##### Threshold grid of biomarker to search and obtain results for
  ####################################################

  # If biomarker is continuous (>10 unique values) then choose grid based on quantiles
  # Otherwise, use all values
  if(length(unique(data_biased[[marker]])) > 15) {
    # Choose grid of 15 thresholds
    # Make sure the thresholds are supported by data (so that A >=v has enough people)
    #  seq(0, 1, length.out = threshold_grid_size) will error code
    # Upper threshold must be less than the 0.99 quantile
    thresh_grid <- unique(quantile(data_biased[[marker]], seq(0, 0.9, length.out = threshold_grid_size), na.rm=T))
  } else {
    thresh_grid <- sort(unique(data_biased[[marker]]))
    #thresh_grid <- thresh_grid[-1]
    print(thresh_grid)
  }

  # Store list mapping nodes in NPSEM to variable names
  node_list <- list(W = baseline, Y = outcome, A  = marker, weights = "wt")
  nodes <- unlist(node_list, use.names = F)


  ####################################################
  # Run thresholdTMLE
  ####################################################


  esttmle <- suppressWarnings(thresholdTMLE(as.data.frame(data_full), node_list, thresholds = thresh_grid,  biased_sampling_strata = "stratum",  biased_sampling_indicator = "sampled", lrnr_A = lrnr, lrnr_Y = lrnr, lrnr_Delta = NULL))
  head(esttmle)
  # Save estimates and CI of threshold-response function
  write.csv(esttmle, file.path(
    folder, paste0("tmleThresh_",   stringr::str_replace(marker, "\\.", "_")
,"_", save_file_key, ".csv")))
  # Clean up table
  esttmle_table <- esttmle
  head(esttmle_table)
  esttmle_table[,1] <- round(esttmle_table[,1], 3)
  esttmle_table[,2] <- round(esttmle_table[,2], 5)
  esttmle_table[,3] <- round(esttmle_table[,3], 5)
  esttmle_table[,4] <- round(esttmle_table[,4], 5)
  esttmle_table[,5] <- round(esttmle_table[,5], 5)
  esttmle_table[,6] <- round(esttmle_table[,6], 5)
  esttmle_table[,7] <- round(esttmle_table[,7], 5)
  esttmle_table <- data.frame(esttmle_table,   paste0(gsub("e\\+0", "$*10^{",format( 10^esttmle_table[,1], scientific = T, digits = 3)) , "}$"))
  esttmle_table[,8] <-  paste0(gsub("e\\-0", "*$10^{-",esttmle_table[,8]), "}$")
  
  
  esttmle_table <- esttmle_table[,c(1,8,2, 4,5, 6,7)]
  esttmle_table[esttmle_table[,3]<0,3] <- 0
  colnames(esttmle_table) <- c("log$_{10}$-Threshold", "Threshold", "Risk estimate", "CI left", "CI right","CI left", "CI right")
  # Save nice latex table
  esttmle_table
  print(esttmle_table)
  index_to_show <-   unique(round(seq.int(1, nrow(esttmle_table), length.out =10)))
  rownames(esttmle_table) <- NULL
  kable(esttmle_table[index_to_show,c(1,2,3,4,5)], row.names = F, format="latex", booktabs=TRUE,linesep = "", escape = F) %>% kable_styling(latex_options = c("scaled_down", "striped"), )%>%
    save_kable(file = file.path(folder,
      paste0("TABLE_",   stringr::str_replace(marker, "\\.", "_")
, "_", save_file_key,"_pointwise.pdf")))
  kable(esttmle_table[index_to_show,c(1,2,3,6,7)], row.names = F, format="latex", booktabs=TRUE,linesep = "", escape = F) %>% kable_styling(latex_options = c("scaled_down", "striped"), )%>%
    save_kable(file = file.path(folder,
      paste0("TABLE_",   stringr::str_replace(marker, "\\.", "_")
,"_", save_file_key, "_simult.pdf")))
  output_list[[marker]] <- esttmle

}
summary_of_results <- output_list

print(summary_of_results)


#######################################################
###### PLOTS of threshold-response and inverse function
#######################################################
library(cowplot)
library(scales)
plot_list <- list()
# For each marker
for(marker in names(output_list)) {
  esttmle <- output_list[[marker]]
  # labels of plot
  main <- "Threshold-response function, "
  labx <- marker

  if(marker == "RSVA.log10cordoverd0") {
    labx <- paste0("RSVA.cordoverd0")
  }
  if(marker == "RSVB.log10cordoverd0") {
    labx <- paste0("RSVB.cordoverd0")
  }
  labx <- gsub("log10", "", marker )
  #######################################################
  # Plot threshold-response with point-wise confidence intervals
  #######################################################
  ident <- function(x) x
  col <- c(col2rgb("olivedrab3")) # orange, darkgoldenrod2
  col <- rgb(col[1], col[2], col[3], alpha=255*0.4, maxColorValue=255)
  #hist(na.omit(data_biased[[marker]]),col=col,axes=F,labels=F,main="",xlab="",ylab="",breaks=10,border=0,freq=F)
  # Get initial threshold-response plot with simultaneous CI
  v <- plot_threshold_response(esttmle, simultaneous_CI = F, monotone = F)
  keep <- data$sampled==1


  data_tmp <- na.omit(data[keep, c(marker, outcome, "wt"), drop = F])
  max_thresh <- max(data_tmp[data_tmp[[outcome]]==1, marker])
  #out <- hist(na.omit(data_tmp[[marker]]), col=col,axes=F,labels=F,main="",xlab="",ylab="",breaks=10,border=0,freq=F)
  scale_coef <- max(v$data$upper) * 1

  #coef <- max(out$density)/scale_coef
  RCDF <- function(a) {
    sum(data_tmp$wt *(data_tmp[[marker]] >=a)) / sum(data_tmp$wt)*scale_coef
  }
  RCDF <- Vectorize(RCDF)
  changeSciNot <- function(n) {
    output <- format(10^n, scientific = TRUE, digits = 2) #Transforms the number into scientific notation even if small
    output <- sub("e", "*10^", output) #Replace e with 10^
    output <- sub("\\+0?", "", output) #Remove + symbol and leading zeros on expoent, if > 1
    output <- sub("-0?", "-", output) #Leaves - symbol but removes leading zeros on expoent, if < 1
    output
  }
  changeSciNot <- Vectorize(changeSciNot)
  
  plot <- v + scale_x_continuous(breaks = unique(c(ceiling(quantile(v$data$cutoffs)), floor(quantile(v$data$cutoffs)) )),
                                 labels = trans_format("ident", math_format(10^.x)), limits =  c(floor(min(esttmle[,1])), ceiling(max(esttmle[,1]))))   + xlab(paste0(marker, " threshold")) +ylab("Probability of RSV disease")+ theme(axis.text.x = element_text(angle = 0, hjust = 1)) + xlab(labx) +
    stat_function(fun = RCDF, color = col,  geom = "area", fill = col, alpha = 0.2) +
    scale_y_continuous(
      name = "Probability of RSV disease",
      sec.axis = sec_axis(~./scale_coef, name="Reverse CDF"), n.breaks = 10
    )+ geom_vline(xintercept=max_thresh, colour="red", linetype = "longdash" )
  #print(ggplot() + geom_histogram(data = data_tmp, aes_string(x = marker, y = "..density../..coef.."), fill = col, colour = col, size =0, alpha = 0.2, bins=10))
  

  #print(ggplot() + geom_histogram(data = data_tmp, aes_string(x = marker, y = "..density../..coef.."), fill = col, colour = col, size =0, alpha = 0.2, bins=10))


  ggsave(filename = file.path(
    folder, paste0("PLOT_",   stringr::str_replace(marker, "\\.", "_")
,"_", save_file_key,"_pointwiseCI", ".pdf")),
    plot = plot, height = 6, width = 7)

  #######################################################
  # Plot threshold-response with point-wise confidence intervals assuming Monotonicity
  #######################################################
  v <- plot_threshold_response(esttmle, simultaneous_CI = F, monotone = T)
keep <- data$sampled==1
  data_tmp <- na.omit(data[keep, c(marker, outcome, "wt"), drop = F])
  max_thresh <- max(data_tmp[data_tmp[[outcome]]==1, marker])
  scale_coef <- max(v$data$upper) * 1
  RCDF <- function(a) {
    sum(data_tmp$wt *(data_tmp[[marker]] >=a)) / sum(data_tmp$wt)*scale_coef
  }
  RCDF <- Vectorize(RCDF)

  changeSciNot <- function(n) {
    output <- format(10^n, scientific = TRUE, digits = 2) #Transforms the number into scientific notation even if small
    output <- sub("e", "*10^", output) #Replace e with 10^
    output <- sub("\\+0?", "", output) #Remove + symbol and leading zeros on expoent, if > 1
    output <- sub("-0?", "-", output) #Leaves - symbol but removes leading zeros on expoent, if < 1
    output
  }
  changeSciNot <- Vectorize(changeSciNot)

  plot <- v + scale_x_continuous(breaks = unique(c(ceiling(quantile(v$data$cutoffs)), floor(quantile(v$data$cutoffs)) )),
                                 labels = trans_format("ident", math_format(10^.x)), limits =  c(floor(min(esttmle[,1])), ceiling(max(esttmle[,1]))))   + xlab(paste0(marker, " threshold")) +ylab("Probability of RSV disease")+ theme(axis.text.x = element_text(angle = 0, hjust = 1)) + xlab(labx) +
    stat_function(fun = RCDF, color = col,  geom = "area", fill = col, alpha = 0.2) +
    scale_y_continuous(
      name = "Probability of RSV disease",
      sec.axis = sec_axis(~./scale_coef, name="Reverse CDF"), n.breaks = 10
    )+ geom_vline(xintercept=max_thresh, colour="red", linetype = "longdash" )
  #print(ggplot() + geom_histogram(data = data_tmp, aes_string(x = marker, y = "..density../..coef.."), fill = col, colour = col, size =0, alpha = 0.2, bins=10))

  ggsave(filename = file.path(
    folder, paste0("PLOT_",   stringr::str_replace(marker, "\\.", "_")
,"_", save_file_key,"_pointwiseCI_monotone", ".pdf")),
    plot = plot, height = 6, width = 7)
  
  

  #######################################################
  # Plot threshold-response with simultaneous confidence intervals
  #######################################################
  head(esttmle)
  v <- plot_threshold_response(esttmle, simultaneous_CI = T)
keep <- data$sampled==1
  data_tmp <- na.omit(data[keep, c(marker, outcome, "wt"), drop = F])
  max_thresh <- max(data_tmp[data_tmp[[outcome]]==1, marker])
  scale_coef <- max(v$data$upper) * 1

  RCDF <- function(a) {
    sum(data_tmp$wt *(data_tmp[[marker]] >=a)) / sum(data_tmp$wt)*scale_coef
  }
  RCDF <- Vectorize(RCDF)

  changeSciNot <- function(n) {
    output <- format(10^n, scientific = TRUE, digits = 2) #Transforms the number into scientific notation even if small
    output <- sub("e", "*10^", output) #Replace e with 10^
    output <- sub("\\+0?", "", output) #Remove + symbol and leading zeros on expoent, if > 1
    output <- sub("-0?", "-", output) #Leaves - symbol but removes leading zeros on expoent, if < 1
    output
  }
  changeSciNot <- Vectorize(changeSciNot)
  
  plot <- v + scale_x_continuous(breaks = unique(c(ceiling(quantile(v$data$cutoffs)), floor(quantile(v$data$cutoffs)) )),
                                 labels = trans_format("ident", math_format(10^.x)), limits =  c(floor(min(esttmle[,1])), ceiling(max(esttmle[,1]))))   + xlab(paste0(marker, " threshold")) +ylab("Probability of RSV disease")+ theme(axis.text.x = element_text(angle = 0, hjust = 1)) + xlab(labx) +
    stat_function(fun = RCDF, color = col,  geom = "area", fill = col, alpha = 0.2) +
    scale_y_continuous(
      name = "Probability of RSV disease",
      sec.axis = sec_axis(~./scale_coef, name="Reverse CDF"), n.breaks = 10
    )+ geom_vline(xintercept=max_thresh, colour="red", linetype = "longdash" )
  #print(ggplot() + geom_histogram(data = data_tmp, aes_string(x = marker, y = "..density../..coef.."), fill = col, colour = col, size =0, alpha = 0.2, bins=10))
  
  ggsave(filename = file.path(
    folder, paste0("PLOT_",   stringr::str_replace(marker, "\\.", "_")
, "_", save_file_key,"_simultaneousCI.pdf")),
    plot = plot, height = 6, width = 7)
  
  
  #######################################################
  # Plot threshold-response with simultaneous confidence intervals assuming monotonicity
  #######################################################
  v <- plot_threshold_response(esttmle, simultaneous_CI = T, monotone = T)
keep <- data$sampled==1
  data_tmp <- na.omit(data[keep, c(marker, outcome, "wt"), drop = F])
  max_thresh <- max(data_tmp[data_tmp[[outcome]]==1, marker])
  scale_coef <- max(v$data$upper) * 1

  RCDF <- function(a) {
    sum(data_tmp$wt *(data_tmp[[marker]] >=a)) / sum(data_tmp$wt)*scale_coef
  }
  RCDF <- Vectorize(RCDF)

  changeSciNot <- function(n) {
    output <- format(10^n, scientific = TRUE, digits = 2) #Transforms the number into scientific notation even if small
    output <- sub("e", "*10^", output) #Replace e with 10^
    output <- sub("\\+0?", "", output) #Remove + symbol and leading zeros on expoent, if > 1
    output <- sub("-0?", "-", output) #Leaves - symbol but removes leading zeros on expoent, if < 1
    output
  }
  changeSciNot <- Vectorize(changeSciNot)
  
  plot <- v + scale_x_continuous(breaks = unique(c(ceiling(quantile(v$data$cutoffs)), floor(quantile(v$data$cutoffs)) )),
                                 labels = trans_format("ident", math_format(10^.x)), limits =  c(floor(min(esttmle[,1])), ceiling(max(esttmle[,1]))))   + xlab(paste0(marker, " threshold")) +ylab("Probability of RSV disease")+ theme(axis.text.x = element_text(angle = 0, hjust = 1)) + xlab(labx) +
    stat_function(fun = RCDF, color = col,  geom = "area", fill = col, alpha = 0.2) +
    scale_y_continuous(
      name = "Probability of RSV disease",
      sec.axis = sec_axis(~./scale_coef, name="Reverse CDF"), n.breaks = 10
    )+ geom_vline(xintercept=max_thresh, colour="red", linetype = "longdash" )
  #print(ggplot() + geom_histogram(data = data_tmp, aes_string(x = marker, y = "..density../..coef.."), fill = col, colour = col, size =0, alpha = 0.2, bins=10))
  
  # v <- plot_threshold_response(esttmle, simultaneous_CI = T)
  # # Save plot of threshold response function
  # plot <- v + scale_x_continuous(n.breaks = 7,
  #   labels = trans_format("ident", math_format(10^.x))) + scale_y_continuous(n.breaks = 10)  + xlab(paste0(marker, " threshold")) +ylab("Probability of RSV disease")+ theme(axis.text.x = element_text(angle = 0, hjust = 1)) + xlab(labx)
  ggsave(filename = file.path(
    folder, paste0("PLOT_",   stringr::str_replace(marker, "\\.", "_")
, "_", save_file_key,"_simultaneousCI_monotone.pdf")),
    plot = plot, height = 6, width = 7)
  
  
#
#   try({
#     #######################################################
#     # Plot INVERSE threshold-response with simultaneous confidence intervals, ASSUMES monotonicity
#     #######################################################
#
#     plot_list[[marker]] <- plot
#
#     # Risks to estimate threshold of protection at.
#     risks <- risks_to_estimate_thresh_of_protection
#     if(is.null(risks)){
#       risks <-  unique(round(seq(max(min(esttmle[esttmle[,2]>0.0005,2]), 0.001), max(esttmle[,2]), length.out = 15),4))
#     }
#     if(length(risks)==1) {
#       risks <- c(risks, risks*2)
#     }
#     plot <- plot_inverse_threshold_response(esttmle, risks = risks)
#     # Save csv file with threshold of protection estimates and CI
#     write.csv(plot$inv_estimates,  file.path(
#       paste0("tmleThresh_inverse_", marker,"_", save_file_key, ".csv")))
#
#     inv_est <- plot$inv_estimates
#     inv_est[,c(2,3,4)] <- round(inv_est[,c(2,3,4)],2)
#     inv_est[,c(1)] <- round(inv_est[,c(1)],4)
#     #inv_est
#
#
#
#     tmp1 <- gsub("e\\+0", "*10^",format( 10^inv_est[,2], scientific = T, digits = 3))
#     tmp2 <- gsub("e\\+0", "*10^",format( 10^inv_est[,4], scientific = T, digits = 3))
#     tmp3 <- gsub("e\\+0", "*10^",format( 10^inv_est[,3], scientific = T, digits = 3))
#     inv_est[,2] <- tmp1
#     inv_est[,4] <- tmp2
#
#     inv_est <- data.frame(inv_est,  tmp3 )
#
#     #inv_est <- data.frame(inv_est, paste0("10^", round(inv_est[,3],1) ))
#     inv_est <- inv_est[,c(1,3, 5,2,4)]
#
#     colnames(inv_est) <- c("Risk", "log10-est.", "Threshold est.", "CI left", "CI right")
#     inv_est <- inv_est[,c(1,3,4,5)]
#
#     # Nice latex table
#     kable(inv_est ,format="latex", booktabs=TRUE) %>% kable_styling(latex_options = c("scaled_down", "striped"))%>%
#       save_kable(file = file.path(
#         paste0("tmleThresh_inverse_TABLE", marker, "_", save_file_key,".pdf")))
#
#     # Save plot of inverse threshold response function
#     plot <- plot$plot + scale_y_continuous(
#       labels = trans_format("ident", math_format(10^.x))) + scale_x_continuous(n.breaks = 7)  + ylab(paste0(marker, " threshold")) +xlab("Probability of RSV disease")+ theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylab(labx) +theme(plot.title = element_text(size=12))
#     ggsave(filename = file.path(
#       paste0("tmleThresh_inverse_", marker, "_", save_file_key,".pdf")),
#       plot = plot, height = 6, width = 7)
#
#
#     #######################################################
#     # Plot INVERSE threshold-response with point-wise confidence intervals, ASSUMES monotonicity
#     ######################################################
#     #### NOTE !!!!!! These results assume monotonicity. This code is the most likely to break if there are too few thresholds estimated or if the function is very weird.
#
#     plot <- plot_inverse_threshold_response(esttmle, risks = risks, simultaneous_CI = F)
#     plot
#     # Save csv file with threshold of protection estimates and CI
#     write.csv(plot$inv_estimates,  file.path(
#       paste0("tmleThresh_inverse_", marker,"_", save_file_key, ".csv")))
#
#     inv_est <- plot$inv_estimates
#     inv_est[,c(2,3,4)] <- round(inv_est[,c(2,3,4)],2)
#     inv_est[,c(1)] <- round(inv_est[,c(1)],4)
#     #inv_est
#
#
#     tmp1 <- gsub("e\\+0", "*10^",format( 10^inv_est[,2], scientific = T, digits = 3))
#     tmp2 <- gsub("e\\+0", "*10^",format( 10^inv_est[,4], scientific = T, digits = 3))
#     tmp3 <- gsub("e\\+0", "*10^",format( 10^inv_est[,3], scientific = T, digits = 3))
#     inv_est[,2] <- tmp1
#     inv_est[,4] <- tmp2
#
#     inv_est <- data.frame(inv_est,  tmp3 )
#
#     #inv_est <- data.frame(inv_est, paste0("10^", round(inv_est[,3],1) ))
#     inv_est <- inv_est[,c(1,3, 5,2,4)]
#
#     colnames(inv_est) <- c("Risk", "log10-est.", "Threshold est.", "CI left", "CI right")
#     inv_est <- inv_est[,c(1,3,4,5)]
#
#     # Nice latex table
#     kable(inv_est ,format="latex", booktabs=TRUE) %>% kable_styling(latex_options = c("scaled_down", "striped"))%>%
#       save_kable(file = file.path(
#         paste0("tmleThresh_inverse_TABLE", marker,"_", save_file_key, "_pointwise.pdf")))
#
#     # Save plot of inverse threshold response function
#     plot <- plot$plot + scale_y_continuous(
#       labels = trans_format("ident", math_format(10^.x))) + scale_x_continuous(n.breaks = 7)  + ylab(paste0(marker, " threshold")) +xlab("Probability of RSV disease")+ theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylab(labx) +theme(plot.title = element_text(size=12))
#     ggsave(filename = file.path(
#       paste0("tmleThresh_inverse_", marker, "_", save_file_key,"_pointwise.pdf")),
#       plot = plot, height = 6, width = 7)
#   })
}

