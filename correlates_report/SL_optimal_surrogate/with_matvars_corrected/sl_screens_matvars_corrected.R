## create SL screens, algorithm + screen combinations

check_include_foldvars <- function(vars, foldvars, vecfoldvars){
  if(length(vecfoldvars)!=0){
    names(vecfoldvars) <- foldvars
    
    if(any(names(vars) %in% "EIA.log10d14overd0") == TRUE){
      if(vars[names(vars) %in% "EIA.log10d14overd0"] == TRUE){
        vecfoldvars[names(vecfoldvars) %in% "EIA.2fold.d14overd0"] <- TRUE
        vecfoldvars[names(vecfoldvars) %in% "EIA.4fold.d14overd0"] <- TRUE
      }
    }
    if(any(names(vars) %in% "PCA.log10d14overd0") == TRUE){
      if(vars[names(vars) %in% "PCA.log10d14overd0"] == TRUE){
        vecfoldvars[names(vecfoldvars) %in% "PCA.2fold.d14overd0"] <- TRUE
        #vecfoldvars[names(vecfoldvars) %in% "PCA.4fold.d14overd0"] <- TRUE
      }
    }
    if(any(names(vars) %in% "RSVA.log10d14overd0") == TRUE){
      if(vars[names(vars) %in% "RSVA.log10d14overd0"] == TRUE){
        vecfoldvars[names(vecfoldvars) %in% "RSVA.2fold.d14overd0"] <- TRUE
        vecfoldvars[names(vecfoldvars) %in% "RSVA.4fold.d14overd0"] <- TRUE
      }
    }
    if(any(names(vars) %in% "RSVB.log10d14overd0") == TRUE){
      if(vars[names(vars) %in% "RSVB.log10d14overd0"] == TRUE){
        vecfoldvars[names(vecfoldvars) %in% "RSVB.2fold.d14overd0"] <- TRUE
        vecfoldvars[names(vecfoldvars) %in% "RSVB.4fold.d14overd0"] <- TRUE
      }
    }
    vars <- c(vars, vecfoldvars)
  }
  return(vars) 
}

## -------------------------------------------------------------------------------------
## SL screens; all models adjust for baseline maternal enrollment variables
## -------------------------------------------------------------------------------------
## screen based on logistic regression univariate p-value < level
rank_univariate_logistic_pval_plus_exposure <- function(Y, X, family, obsWeights, id, ...) {
  ## logistic regression of outcome on each variable
  listp <- apply(X, 2, function(x, Y, family) {
    summ <- coef(summary(glm(Y ~ x + X$age.at.trt + X$age.at.trt.cat + X$bmi + X$mhsmopr + X$m.ast + X$child5 + X$season + X$smoker + X$daycare, 
                             family = family, weights = obsWeights)))
    ifelse(dim(summ)[1] > 1, summ[2, 4], 1)
  }, Y = Y, family = family)
  ## rank the p-values; give all maternal enrollment variables the lowest rank (will always set to TRUE anyways)
  ranked_vars <- rank(listp, ties = "average")
  ranked_vars[names(X) %in% c("age.at.trt", "age.at.trt.cat", "bmi", "mhsmopr", "m.ast", "child5", "season", "smoker", "daycare")] <- 999
  return(ranked_vars)
}

# no screen (only the top 4 markers with lowest univariate logistic pvalue + maternal enrollment variables)
screen_all_plus_exposure<- function(Y, X, family, obsWeights, id, nVar=6,...){
  # implement screens upon excluding the indicator fold variables
  foldvars <- names(X[names(X) %in% c("EIA.2fold.d14overd0", "EIA.4fold.d14overd0", "PCA.2fold.d14overd0", "PCA.4fold.d14overd0", "RSVA.2fold.d14overd0",
                                      "RSVA.4fold.d14overd0", "RSVB.2fold.d14overd0", "RSVB.4fold.d14overd0")])
  X_nofoldvars <- X[!names(X) %in% foldvars]
  
  #X contain baseline variables
  ## logistic regression of outcome on each variable
  vars <- rep(TRUE, ncol(X_nofoldvars))
  listp <- apply(X_nofoldvars, 2, function(x, Y, family) {
    summ <- coef(summary(glm(Y ~ x + X$age.at.trt + X$age.at.trt.cat + X$bmi + X$mhsmopr + X$m.ast + X$child5 + X$season + X$smoker + X$daycare, 
                             family = family, weights = obsWeights)))
    ifelse(dim(summ)[1] > 1, summ[2, 4], 1)
  }, Y = Y, family = family)
  
  nomatvars <- listp[!names(listp) %in% c("age.at.trt", "age.at.trt.cat", "bmi", "mhsmopr", "m.ast", "child5", "season", "smoker", "daycare")]
  ranked_nomatvars <- rank(nomatvars, ties = "average")
  nomatvars <- rep(TRUE, length(nomatvars))
  nomatvars[nomatvars][ranked_nomatvars > nVar] <- FALSE
  matvars <- rep(TRUE, length(c("age.at.trt", "age.at.trt.cat", "bmi", "mhsmopr", "m.ast", "child5", "season", "smoker", "daycare"))) 
  vars <- c(matvars, nomatvars)
  names(vars) <- names(X_nofoldvars)
  
  # include the indicator fold variables. Make them TRUE only if corresponding quantitative fold variable is TRUE
  vecfoldvars <- rep(FALSE, ncol(X[names(X) %in% foldvars]))
  vars <- check_include_foldvars(vars, foldvars, vecfoldvars)
  return(vars)
}

## screen based on lasso 
screen_glmnet_plus_exposure <- function(Y, X, family, obsWeights, id, alpha = 1, minscreen = 2, nfolds = 10, nlambda = 100, nVar = 6, ...) {
  # implement screens upon excluding the indicator fold variables
  foldvars <- names(X[names(X) %in% c("EIA.2fold.d14overd0", "EIA.4fold.d14overd0", "PCA.2fold.d14overd0", "PCA.4fold.d14overd0", "RSVA.2fold.d14overd0",
                                      "RSVA.4fold.d14overd0", "RSVB.2fold.d14overd0", "RSVB.4fold.d14overd0")])
  X_nofoldvars <- X[!names(X) %in% foldvars]
  
  vars <- screen.glmnet(Y, X_nofoldvars, family, obsWeights, id, alpha = 1, minscreen = 2, nfolds = 10, nlambda = 100, ...)
  # also keep the first nine columns of X (correspond to maternal enrollment variables)
  vars[names(X_nofoldvars) %in% c("age.at.trt", "age.at.trt.cat", "bmi", "mhsmopr", "m.ast", "child5", "season", "smoker", "daycare")] <- TRUE
  # keep only a max of nVar immune markers; rank by univariate p-value in a model adjusting for maternal enrollment variables
  X_initial_screen <- X_nofoldvars %>%
    select(names(X_nofoldvars)[vars], "age.at.trt", "age.at.trt.cat", "bmi", "mhsmopr", "m.ast", "child5", "season", "smoker", "daycare")
  ranked_vars <- rank_univariate_logistic_pval_plus_exposure(Y, X_initial_screen, family, obsWeights, id)
  
  # re-rank without including maternal enrollment variables
  ranked_vars <- c(ranked_vars[names(ranked_vars) %in% c("age.at.trt", "age.at.trt.cat", "bmi", "mhsmopr", "m.ast", "child5", "season", "smoker", "daycare")],
                   rank(ranked_vars[!names(ranked_vars) %in% c("age.at.trt", "age.at.trt.cat", "bmi", "mhsmopr", "m.ast", "child5", "season", "smoker", "daycare")]))
  
  vars[vars][ranked_vars > nVar] <- FALSE
  # also keep the first nine columns of X (correspond to maternal enrollment variables)
  vars[names(X_nofoldvars) %in% c("age.at.trt", "age.at.trt.cat", "bmi", "mhsmopr", "m.ast", "child5", "season", "smoker", "daycare")] <- TRUE
  
  # include the indicator fold variables. Make them TRUE only if corresponding quantitative fold variable is TRUE
  vecfoldvars <- rep(FALSE, ncol(X[names(X) %in% foldvars]))
  vars <- check_include_foldvars(vars, foldvars, vecfoldvars)
  return(vars)
}

## screen based on logistic regression univariate p-value < level
screen_univariate_logistic_pval_plus_exposure <- function(Y, X, family, obsWeights, id, minPvalue = 0.1, minscreen = 2, nVar = 6, ...) {
  # implement screens upon excluding the indicator fold variables
  foldvars <- names(X[names(X) %in% c("EIA.2fold.d14overd0", "EIA.4fold.d14overd0", "PCA.2fold.d14overd0", "PCA.4fold.d14overd0", "RSVA.2fold.d14overd0",
                                      "RSVA.4fold.d14overd0", "RSVB.2fold.d14overd0", "RSVB.4fold.d14overd0")])
  X_nofoldvars <- X[!names(X) %in% foldvars]
  
  ## logistic regression of outcome on each variable
  listp <- apply(X_nofoldvars, 2, function(x, Y, family) {
    summ <- coef(summary(glm(Y ~ x + X$age.at.trt + X$age.at.trt.cat + X$bmi + X$mhsmopr + X$m.ast + X$child5 + X$season + X$smoker + X$daycare, 
                             family = family, weights = obsWeights)))
    ifelse(dim(summ)[1] > 1, summ[2, 4], 1)
  }, Y = Y, family = family)
  vars <- (listp <= minPvalue)
  # also keep the first nine columns of X (correspond to maternal enrollment variables)
  vars[names(X_nofoldvars) %in% c("age.at.trt", "age.at.trt.cat", "bmi", "mhsmopr", "m.ast", "child5", "season", "smoker", "daycare")] <- TRUE
  if (sum(vars) < minscreen) {
    warning("number of variables with p value less than minPvalue is less than minscreen")
    vars[rank(listp) <= minscreen] <- TRUE
  }
  # keep only a max of nVar immune markers; rank by univariate p-value in a model adjusting for maternal enrollment variables
  X_initial_screen <- X_nofoldvars %>%
    select(names(X_nofoldvars)[vars], "age.at.trt", "age.at.trt.cat", "bmi", "mhsmopr", "m.ast", "child5", "season", "smoker", "daycare")
  ranked_vars <- rank_univariate_logistic_pval_plus_exposure(Y, X_initial_screen, family, obsWeights, id)
  
  # re-rank without including maternal enrollment variables
  ranked_vars <- c(ranked_vars[names(ranked_vars) %in% c("age.at.trt", "age.at.trt.cat", "bmi", "mhsmopr", "m.ast", "child5", "season", "smoker", "daycare")],
                   rank(ranked_vars[!names(ranked_vars) %in% c("age.at.trt", "age.at.trt.cat", "bmi", "mhsmopr", "m.ast", "child5", "season", "smoker", "daycare")]))
  
  vars[vars][ranked_vars > nVar] <- FALSE
  # also keep the first nine columns of X (correspond to maternal enrollment variables)
  vars[names(X_nofoldvars) %in% c("age.at.trt", "age.at.trt.cat", "bmi", "mhsmopr", "m.ast", "child5", "season", "smoker", "daycare")] <- TRUE
  
  # include the indicator fold variables. Make them TRUE only if corresponding quantitative fold variable is TRUE
  vecfoldvars <- rep(FALSE, ncol(X[names(X) %in% foldvars]))
  vars <- check_include_foldvars(vars, foldvars, vecfoldvars)
  return(vars)
}


screen_highcor_random_plus_exposure <- function(Y, X, family, obsWeights, id, nVar = 6, ...) {
  
  # implement screens upon excluding the indicator fold variables
  foldvars <- names(X[names(X) %in% c("EIA.2fold.d14overd0", "EIA.4fold.d14overd0", "PCA.2fold.d14overd0", "PCA.4fold.d14overd0", "RSVA.2fold.d14overd0",
                                      "RSVA.4fold.d14overd0", "RSVB.2fold.d14overd0", "RSVB.4fold.d14overd0")])
  X_nofoldvars <- X[!names(X) %in% foldvars]
  
  # set all vars to FALSE
  vars <- rep(FALSE, ncol(X_nofoldvars))
  # compute pairwise correlations between all marker vars
  cors <- cor(X_nofoldvars, method = "spearman")
  diag(cors) <- NA
  cor_less_0.9 <- (cors <= 0.9)
  # screen out those with r > 0.9
  vars <- apply(cor_less_0.9, 1, function(x) all(x, na.rm = TRUE))
  
  # test with cor less than 0.2 !!!
  # cor_less_0.2 <- (cors <= 0.2)
  # vars <- apply(cor_less_0.2, 1, function(x) all(x, na.rm = TRUE))
  
  # if cor is greater than 0.9 for any pair of variables, pick up one of the variables at random!
  cormat = cor_less_0.9
  long.cormat = data.frame(row=rownames(cormat)[row(cormat)[upper.tri(cormat)]],       # gets only upper triangle of symmetrical corr matrix and puts data in long format
                           col=colnames(cormat)[col(cormat)[upper.tri(cormat)]], 
                           corr=cormat[upper.tri(cormat)]) %>%
    filter(corr == "FALSE") 
  
  if(dim(long.cormat)[1] > 0) {
    # select random element out of any pair
    long.cormat$randCol = apply(long.cormat, 1, function(x)  sample(c(x[1], x[2]), 1, replace=T))
    # get the unique columns
    randCols = unique(long.cormat$randCol)
    # # just select 2 for now for checking!
    # randCols = randCols[1:2]
    vars[randCols] <- TRUE
  }
  
  # keep only a max of nVar immune markers; select randomly
  nomatvars <- vars[!names(vars) %in% c("age.at.trt", "age.at.trt.cat", "bmi", "mhsmopr", "m.ast", "child5", "season", "smoker", "daycare")]
  if(length(nomatvars[nomatvars]) > 4){
    sampCols <- sample(names(nomatvars[nomatvars]), nVar)
  }else{
    sampCols <- names(nomatvars[nomatvars])
  }
  
  vars[!names(vars) %in% sampCols] <- FALSE
  
  ## make sure that maternal enrollment vars are true
  vars[names(X_nofoldvars) %in% c("age.at.trt", "age.at.trt.cat", "bmi", "mhsmopr", "m.ast", "child5", "season", "smoker", "daycare")] <- TRUE
  
  # include the indicator fold variables. Make them TRUE only if corresponding quantitative fold variable is TRUE
  vecfoldvars <- rep(FALSE, ncol(X[names(X) %in% foldvars]))
  vars <- check_include_foldvars(vars, foldvars, vecfoldvars)
  
  return(vars)
}


# screens <- c("screen_all_plus_exposure","screen_glmnet_plus_exposure", 
#              "screen_univariate_logistic_pval_plus_exposure",
#              "screen_highcor_random_plus_exposure")


# ## -------------------------------------------------------------------------------------
# ## assign the screen function with different maximum covariates to the envirmment, k from 1 to 8
# ## -------------------------------------------------------------------------------------
# 
# lapply(1:8,function(i){
#   assign(x=paste0('screen_all_plus_exposure_',i),value=
#            eval(parse(text=paste0('function(Y, X, family, obsWeights, id, nVar=',i,',...){screen_all_plus_exposure(Y, X, family, obsWeights, id, nVar=',i, ',...)}'))),envir = .GlobalEnv) }
# )
# 
# lapply(1:8,function(i){
#   assign(x=paste0('screen_glmnet_plus_exposure_',i),value=
#            eval(parse(text=paste0('function(Y, X, family, obsWeights, id, alpha = 1, minscreen = 2, nfolds = 10, nlambda = 100, nVar =',i, ',...){screen_glmnet_plus_exposure(Y, X, family, obsWeights, id, alpha = 1,  nfolds = 10, nlambda = 100, nVar =',i, ',...)}'))),envir = .GlobalEnv) }
# )
# 
# 
# lapply(1:8,function(i){
#   assign(x=paste0('screen_univariate_logistic_pval_plus_exposure_',i),value=
#            eval(parse(text=paste0('function(Y, X, family, obsWeights, id, alpha = 1, minscreen = 2, nfolds = 10, nlambda = 100, nVar =',i, ',...){screen_univariate_logistic_pval_plus_exposure (Y, X, family, obsWeights, id, minPvalue = 0.1, nVar =',i, ',...)}'))),envir = .GlobalEnv) }
# )
# 
# lapply(1:8,function(i){
#   assign(x=paste0('screen_highcor_random_plus_exposure_',i),value=
#            eval(parse(text=paste0('function(Y, X, family, obsWeights, id, alpha = 1, minscreen = 2, nfolds = 10, nlambda = 100, nVar =',i, ',...){screen_highcor_random_plus_exposure(Y, X, family, obsWeights, id, nVar =',i,', ...)}'))),envir = .GlobalEnv) }
# )
## -------------------------------------------------------------------------------------
## SL algorithms
## -------------------------------------------------------------------------------------

## --------------------------------------------------------------------------
## define wrappers that are less memory-intensive than the usual SL functions
## --------------------------------------------------------------------------
# skinny glm
SL.glm.skinny <- function(Y, X, newX, family, obsWeights, ...){
  SL.glm.fit <- SL.glm(Y = Y, X = X, newX = newX, family = family, obsWeights = obsWeights, ...)
  SL.glm.fit$fit$object$y <- NULL
  SL.glm.fit$fit$object$model <- NULL
  SL.glm.fit$fit$object$residuals <- NULL
  SL.glm.fit$fit$object$fitted.values <- NULL
  SL.glm.fit$fit$object$effects <- NULL
  SL.glm.fit$fit$object$qr$qr <- NULL
  SL.glm.fit$fit$object$linear.predictors <- NULL
  SL.glm.fit$fit$object$weights <- NULL
  SL.glm.fit$fit$object$prior.weights <- NULL
  SL.glm.fit$fit$object$data <- NULL
  SL.glm.fit$fit$object$family$variance <- NULL
  SL.glm.fit$fit$object$family$dev.resids <- NULL
  SL.glm.fit$fit$object$family$aic <- NULL
  SL.glm.fit$fit$object$family$validmu <- NULL
  SL.glm.fit$fit$object$family$simulate <- NULL
  attr(SL.glm.fit$fit$object$terms, ".Environment") <- NULL
  attr(SL.glm.fit$fit$object$formula, ".Environment") <- NULL
  return(SL.glm.fit)
}

## skinny glm with interactions
SL.glm.interaction.skinny <- function(Y, X, newX, family, obsWeights, ...){
  SL.glm.fit <- SL.glm.interaction(Y = Y, X = X, newX = newX, family = family, obsWeights = obsWeights, ...)
  SL.glm.fit$fit$object$y <- NULL
  SL.glm.fit$fit$object$model <- NULL
  SL.glm.fit$fit$object$residuals <- NULL
  SL.glm.fit$fit$object$fitted.values <- NULL
  SL.glm.fit$fit$object$effects <- NULL
  SL.glm.fit$fit$object$qr$qr <- NULL
  SL.glm.fit$fit$object$linear.predictors <- NULL
  SL.glm.fit$fit$object$weights <- NULL
  SL.glm.fit$fit$object$prior.weights <- NULL
  SL.glm.fit$fit$object$data <- NULL
  SL.glm.fit$fit$object$family$variance <- NULL
  SL.glm.fit$fit$object$family$dev.resids <- NULL
  SL.glm.fit$fit$object$family$aic <- NULL
  SL.glm.fit$fit$object$family$validmu <- NULL
  SL.glm.fit$fit$object$family$simulate <- NULL
  attr(SL.glm.fit$fit$object$terms, ".Environment") <- NULL
  attr(SL.glm.fit$fit$object$formula, ".Environment") <- NULL
  return(SL.glm.fit)
}

# skinny stepwise with interactions
SL.step.interaction.skinny <- function(Y, X, newX, family, obsWeights, ...){
  SL.step.interaction.fit <- SL.step.interaction(Y = Y, X = X, newX = newX, family = family,
                                                 obsWeights = obsWeights, direction = "forward", ...)
  SL.step.interaction.fit$fit$object$y <- NULL
  SL.step.interaction.fit$fit$object$model <- NULL
  SL.step.interaction.fit$fit$object$residuals <- NULL
  SL.step.interaction.fit$fit$object$fitted.values <- NULL
  SL.step.interaction.fit$fit$object$effects <- NULL
  SL.step.interaction.fit$fit$object$qr$qr <- NULL
  SL.step.interaction.fit$fit$object$linear.predictors <- NULL
  SL.step.interaction.fit$fit$object$weights <- NULL
  SL.step.interaction.fit$fit$object$prior.weights <- NULL
  SL.step.interaction.fit$fit$object$data <- NULL
  SL.step.interaction.fit$fit$object$family$variance <- NULL
  SL.step.interaction.fit$fit$object$family$dev.resids <- NULL
  SL.step.interaction.fit$fit$object$family$aic <- NULL
  SL.step.interaction.fit$fit$object$family$validmu <- NULL
  SL.step.interaction.fit$fit$object$family$simulate <- NULL
  attr(SL.step.interaction.fit$fit$object$terms, ".Environment") <- NULL
  attr(SL.step.interaction.fit$fit$object$formula, ".Environment") <- NULL
  return(SL.step.interaction.fit)
}

# skinny stepwise (forward)
SL.step.skinny <- function(Y, X, newX, family, obsWeights, ...){
  SL.step.fit <- SL.step(Y = Y, X = X, newX = newX, family = family, obsWeights = obsWeights, direction = "forward", ...)
  SL.step.fit$fit$object$y <- NULL
  SL.step.fit$fit$object$model <- NULL
  SL.step.fit$fit$object$residuals <- NULL
  SL.step.fit$fit$object$fitted.values <- NULL
  SL.step.fit$fit$object$effects <- NULL
  SL.step.fit$fit$object$qr$qr <- NULL
  SL.step.fit$fit$object$linear.predictors <- NULL
  SL.step.fit$fit$object$weights <- NULL
  SL.step.fit$fit$object$prior.weights <- NULL
  SL.step.fit$fit$object$data <- NULL
  SL.step.fit$fit$object$family$variance <- NULL
  SL.step.fit$fit$object$family$dev.resids <- NULL
  SL.step.fit$fit$object$family$aic <- NULL
  SL.step.fit$fit$object$family$validmu <- NULL
  SL.step.fit$fit$object$family$simulate <- NULL
  attr(SL.step.fit$fit$object$terms, ".Environment") <- NULL
  attr(SL.step.fit$fit$object$formula, ".Environment") <- NULL
  return(SL.step.fit)
}

# boosted decision stumps
SL.stumpboost <- function(Y, X, newX, family, obsWeights, ...){
  fit <- SL.xgboost(Y = Y, X = X, newX = newX, family = family, obsWeights = obsWeights,
                    max_depth = 1, # so it's only a stump
                    ...)
  return(fit)
}

# naive bayes wrapper
SL.naivebayes <- function(Y, X, newX, family, obsWeights, laplace = 0, ...){
  SuperLearner:::.SL.require("e1071")
  if(family$family == "gaussian"){
    stop("SL.naivebayes only works with binary outcomes")
  }else{
    nb <- naiveBayes(y = Y, x = X, laplace = laplace)
    pred <- predict(nb, newX, type = "raw")[,2]
    out <- list(fit = list(object = nb), pred = pred)
    class(out$fit) <- "SL.naivebayes"
    return(out)
  }
}

# predict method for naive bayes wrapper
predict.SL.naivebayes <- function(object, newdata, ...){
  pred <- predict(object$object, newdata = newdata, type = "raw")[,2]
  return(pred)
}


## -------------------------------------------------------------------------------------
## Add all alg/screen combinations to global environment, create SL library
## -------------------------------------------------------------------------------------

#' This function takes a super learner method wrapper and a super learner
#' screen wrapper and combines them into a single wrapper and makes that
#' wrapper available in the specified environment. It also makes a predict
#' method available in the specified environment.
#' @param method A super learner method wrapper. See ?SuperLearner::listWrappers(what = "method").
#' @param screen A super learner method wrapper. See ?SuperLearner::listWrappers(what = "screen").
#' @param envir The environment to assign the functions to (default is global environment)
#' @param verbose Print a message with the function names confirming their assignment?
assign_combined_function <- function(method, screen, envir = .GlobalEnv,
                                     verbose = TRUE){
  fn <- eval(parse(text =
                     paste0("function(Y, X, newX, obsWeights, family, ...){ \n",
                            "screen_call <- ", screen, "(Y = Y, X = X, newX = newX, obsWeights = obsWeights, family = family, ...) \n",
                            "method_call <- ", method, "(Y = Y, X = X[,screen_call,drop=FALSE], newX = newX[,screen_call,drop = FALSE], obsWeights = obsWeights, family = family, ...) \n",
                            "pred <- method_call$pred \n",
                            "fit <- list(object = method_call$fit$object, which_vars = screen_call) \n",
                            "class(fit) <- paste0('", screen, "', '_', '", method, "') \n",
                            "out <- list(fit = fit, pred = pred) \n",
                            "return(out) \n",
                            "}")))
  fn_name <- paste0(screen,"_",method)
  assign(x = fn_name, value = fn, envir = envir)
  if(verbose){
    message(paste0("Function ", fn_name, " now available in requested environment."))
  }
  if (method == "SL.glmnet") {
    pred_fn <- eval(parse(text =
                            paste0("function(object, newdata, ...){ \n",
                                   "screen_newdata <- newdata[,object$which_vars,drop = FALSE] \n",
                                   "pred <- predict(object$object, type = 'response', newx = as.matrix(screen_newdata), s = 'lambda.min', ...) \n",
                                   "return(pred) \n",
                                   "}")))
  } else if (method == "SL.stumpboost") {
    pred_fn <- eval(parse(text =
                            paste0("function(object, newdata, ...){ \n",
                                   "screen_newdata <- newdata[,object$which_vars,drop = FALSE] \n",
                                   "screen_newdata_2 <- matrix(unlist(lapply(screen_newdata, as.numeric)), nrow=nrow(screen_newdata), ncol=ncol(screen_newdata)) \n",
                                   "pred <- predict(object$object, newdata = screen_newdata_2, ...) \n",
                                   "return(pred) \n",
                                   "}")))
  } else if (method == "SL.naivebayes") {
    pred_fn <- eval(parse(text =
                            paste0("function(object, newdata, ...){ \n",
                                   "screen_newdata <- newdata[,object$which_vars,drop = FALSE] \n",
                                   'pred <- predict(object$object, newdata = screen_newdata, type = "raw", ...)[,2] \n',
                                   "return(pred) \n",
                                   "}")))
  } else if (method == "SL.randomForest") {
    pred_fn <- eval(parse(text =
                            paste0("function(object, newdata, ...){ \n",
                                   "screen_newdata <- newdata[,object$which_vars,drop = FALSE] \n",
                                   "if (object$object$type != 'classification') {
                                    pred <- predict(object$object, newdata = screen_newdata, type = 'response')
                                }else {
                                    pred <- predict(object$object, newdata = screen_newdata, type = 'vote')[,
                                        2]
                                }
                                pred",
                                   "}")))
  }
  else if(method=='SL.mean'){
    pred_fn <- eval(parse(text =
                            paste0("function(object, newdata, ...){ \n",
                                   'pred <- rep.int(object$object, times = nrow(newdata)) \n',
                                   'pred}')))
  }
  else if(method=='SL.cforest'){
    pred_fn<-eval(parse(text =
                          paste0("function(object, newdata, ...){ \n",
                                 'pred <- predict(object=object$object, newdata=newdata[,object$which_vars,drop = FALSE]) \n',
                                 'pred}')))
  }
  else if(method=='SL.xgboost'){
    pred_fn<-eval(parse(text =
                          paste0("function(object, newdata, ...){ \n",
                                 "if (!is.matrix(newdata)) { \n",
                                   "newdata = model.matrix(~. - 1, newdata)}\n",
                                 'pred <- predict(object=object$object, newdata=newdata[,object$which_vars,drop = FALSE]) \n',
                                 'pred}')))
  }
  
  else {
    pred_fn <- eval(parse(text =
                            paste0("function(object, newdata, ...){ \n",
                                   "screen_newdata <- newdata[,object$which_vars,drop = FALSE] \n",
                                   "pred <- predict(object$object, type = 'response', newdata = screen_newdata, ...) \n",
                                   "return(pred) \n",
                                   "}")))
  }
  
  pred_fn_name <- paste0("predict.",screen,"_",method)
  assign(x = pred_fn_name, value = pred_fn, envir = envir)
  if(verbose){
    message(paste0("Function ", pred_fn_name, " now available in requested environment."))
  }
}


# ## -------------------------------------------------------------------------------------
# ## This part is risk score. Add all alg/screen combinations to global environment, create SL library
# ## -------------------------------------------------------------------------------------
# construct_screen <- function(x){
#   #learners in the method1 are also combined with no screen 
#   methods1 <- c("SL.mean","SL.glm","SL.glm.interaction","SL.bayesglm")
#   
#   #learnes in the method2 are learners can have creens
#   methods2 <- c("SL.glm","SL.glm.interaction","SL.bayesglm", "SL.step", "SL.glmnet","SL.gam","SL.cforest","SL.xgboost")
#   
#   screens1 = paste0('screen_all_plus_exposure_',x)
#   screens2 <- paste0(c("screen_glmnet_plus_exposure_", 
#                        "screen_univariate_logistic_pval_plus_exposure_",
#                        "screen_highcor_random_plus_exposure_"),x)
#   
#   screen_method_frame1 <- expand.grid(screen = screens1, method = methods1)
#   screen_method_frame2 <- expand.grid(screen = screens2, method = methods2)
#   
#   apply(screen_method_frame1, 1, function(x) {assign_combined_function(screen = x[1], method = x[2], verbose = FALSE)})
#   apply(screen_method_frame2, 1, function(x) {assign_combined_function(screen = x[1], method = x[2], verbose = FALSE)})
#   
#   SL_library1 <- c(apply(screen_method_frame1, 1, paste0, collapse = "_"))
#   SL_library2 <- c(apply(screen_method_frame2, 1, paste0, collapse = "_"))
#   SL_library <- c(SL_library1, SL_library2)
# }
# 
# SL_library_all = lapply(1:8, construct_screen)
# lapply(1:8,function(x) assign(x=paste0('SL_library_',x),value=SL_library_all[[x]],envir = .GlobalEnv))
# SL_library_1=SL_library_1[c(1:16,20:28)]
# 
# 

## -------------------------------------------------------------------------------------
## This part is objective 3. Add all alg/screen combinations to global environment, create SL library
## -------------------------------------------------------------------------------------

#learners in the method1 are also combined with no screen 
methods1 <- c("SL.mean","SL.glm","SL.glm.interaction","SL.bayesglm")
  
#learnes in the method2 are learners can have creens
methods2 <- c("SL.glm","SL.glm.interaction","SL.bayesglm", "SL.step", "SL.glmnet","SL.gam","SL.cforest","SL.xgboost")
  
screens1 = "screen_all_plus_exposure"
screens2 <- c("screen_glmnet_plus_exposure", 
              "screen_univariate_logistic_pval_plus_exposure",
              "screen_highcor_random_plus_exposure")
  
screen_method_frame1 <- expand.grid(screen = screens1, method = methods1)
screen_method_frame2 <- expand.grid(screen = screens2, method = methods2)
  
apply(screen_method_frame1, 1, function(x) {assign_combined_function(screen = x[1], method = x[2], verbose = FALSE)})
apply(screen_method_frame2, 1, function(x) {assign_combined_function(screen = x[1], method = x[2], verbose = FALSE)})
  
SL_library1 <- c(apply(screen_method_frame1, 1, paste0, collapse = "_"))
SL_library2 <- c(apply(screen_method_frame2, 1, paste0, collapse = "_"))
SL_library <- c(SL_library1, SL_library2)




