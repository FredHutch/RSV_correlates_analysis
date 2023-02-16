## ----package-loading-and-options, warning=FALSE, include=FALSE--------------------------------------------------------------------------------------------------------

## ----load-all-SLobjects, message=FALSE, error=FALSE, warning=FALSE----------------------------------------------------------------------------------------------------
library("cvAUC")
library("conflicted")
library("tidyverse")
library("dplyr")
library("cowplot")
library("ggplot2")
library("vimp")
library("kyotil")
library(gridExtra)
library(cowplot)
library(RSVcorr)
library(here)
source(paste0(here(), "/correlates_report/objective3/with_allvars_n6/define-screens-and-algs.R"))
source(paste0(here(), "/correlates_report/objective3/with_allvars_n6/utils.R"))
#source("plot_assays.R")
source(paste0(here(), "/correlates_report/objective3/with_allvars_n6/make_forest_plot.R"))
method <- "method.CC_nloglik" # since SuperLearner relies on this to be in GlobalEnv
ggplot2::theme_set(theme_cowplot())

load(file = "N:/cavd/Objective 4/GH-VAP/ID33-Wilson/Data/adata/objective3/results_with_allvars_n6/y1_vaccine/y1_vaccine.rda")
load(file = "N:/cavd/Objective 4/GH-VAP/ID33-Wilson/Data/adata/objective3/results_with_allvars_n6/y1_placebo/y1_placebo.rda")
load(file = "N:/cavd/Objective 4/GH-VAP/ID33-Wilson/Data/adata/objective3/results_with_allvars_n6/y2_vaccine/y2_vaccine.rda")
load(file = "N:/cavd/Objective 4/GH-VAP/ID33-Wilson/Data/adata/objective3/results_with_allvars_n6/y2_placebo/y2_placebo.rda")

# filter() is used in both dplyr and stats, so need to set the preference to dplyr
conflict_prefer("filter", "dplyr")
conflict_prefer("summarise", "dplyr")


## ----learner-screens, warning=kable_warnings--------------------------------------------------------------------------------------------------------------------------
caption <- "All learner-screen combinations (28 in total) used as input to the superlearner."

tab <- y2_placebo %>%
  filter(!Learner %in% c("SL", "Discrete SL")) %>%
  filter(!file %in% c("1_only_matvars")) %>%
  select(Learner, Screen) %>%
  mutate(Screen = fct_relevel(Screen, c("all", "glmnet", "univar_logistic_pval",
                                        "highcor_random")),
         Learner = as.factor(Learner),
         Learner = fct_relevel(Learner, c("SL.mean", "SL.glm", "SL.glm.interaction", 
                                          "SL.bayesglm", "SL.step", "SL.glmnet",
                                          "SL.gam", "SL.cforest", "SL.xgboost"))) %>%
  arrange(Learner, Screen) %>% 
  distinct(Learner, Screen) %>%
  rename("Screen*" = Screen) 

tab %>% write.csv(paste0(here(), "/correlates_report/objective3/input/learner-screens.csv"))

## ----All 29 variable sets --------------------------------------------------------------------------------------------------------------------
caption <- "The 29 variable sets on which an estimated optimal surrogate was built."

tab <- data.frame(`Variable Set Name` = c("1_only_matvars",  
                                          "2_varset_EIA.log10d0", "3_varset_EIA.log10d14", "4_varset_EIA.log10d14overd0", "5_varset_EIA.log10cord",
                                          "6_varset_PCA.log10d0", "7_varset_PCA.log10d14", "8_varset_PCA.log10d14overd0", "9_varset_PCA.log10cord", 
                                          "10_varset_RSVA.log10d0", "11_varset_RSVA.log10d14", "12_varset_RSVA.log10d14overd0", "13_varset_RSVA.log10cord",
                                          "14_varset_RSVB.log10d0", "15_varset_RSVB.log10d14", "16_varset_RSVB.log10d14overd0", "17_varset_RSVB.log10cord", 
                                          "18_varset_d0_all4","19_varset_d14_all4", "20_varset_d14overd0_all4", "21_varset_cord_all4",
                                          "22_varset_EIA_all4", "23_varset_PCA_all4", "24_varset_RSVA_all4", "25_varset_RSVB_all4",
                                          "26_varset_d0_d14", "27_varset_d0_foldrise", "28_varset_d0_cord", "29_varset_all16"), 
                  `Variables included in the set` = c("Maternal baseline demographic covariates only (Reference model)",  
                                                      "Maternal covariates + EIA marker at Day 0", 
                                                      "Maternal covariates + EIA marker at Day 14", 
                                                      "Maternal covariates + EIA marker fold-rise", 
                                                      "Maternal covariates + EIA marker at birth (cord)",
                                                      "Maternal covariates + PCA marker at Day 0", 
                                                      "Maternal covariates + PCA marker at Day 14", 
                                                      "Maternal covariates + PCA marker fold-rise", 
                                                      "Maternal covariates + PCA marker at birth (cord)",
                                                      "Maternal covariates + RSVA marker at Day 0", 
                                                      "Maternal covariates + RSVA marker at Day 14", 
                                                      "Maternal covariates + RSVA marker fold-rise", 
                                                      "Maternal covariates + RSVA marker at birth (cord)",
                                                      "Maternal covariates + RSVB marker at Day 0", 
                                                      "Maternal covariates + RSVB marker at Day 14", 
                                                      "Maternal covariates + RSVB marker fold-rise", 
                                                      "Maternal covariates + RSVB marker at birth (cord)",
                                                      "Maternal covariates + all 4 markers at Day 0",
                                                      "Maternal covariates + all 4 markers at Day 14", 
                                                      "Maternal covariates + all 4 markers fold-rise", 
                                                      "Maternal covariates + all 4 markers at birth (cord)",
                                                      "Maternal covariates + all 4 EIA markers (Days 0, 14, fold-rise and cord)", 
                                                      "Maternal covariates + all 4 PCA markers", 
                                                      "Maternal covariates + all 4 RSVA markers", 
                                                      "Maternal covariates + all 4 RSVB markers",
                                                      "Maternal covariates + all 4 markers at Day 0 and Day 14", 
                                                      "Maternal covariates + all 4 markers at Day 0 and their fold-rise", 
                                                      "Maternal covariates + all 4 markers at Day 0 and birth (cord)", 
                                                      "Maternal covariates + all 4 markers at Day 0, 14, at birth and their fold-rise"))

tab %>% write.csv(paste0(here(), "/correlates_report/objective3/input/varsets.csv"))


## ----SLperformance-vacc-y1, warning=kable_warnings--------------------------------------------------------------------------------------------------------------------
caption <- "Performance of Superlearner and the top-performing learner-screen combinations (CV-AUCs with 95\\% CIs) for each of the 29 variable sets using vaccine group and endpoint 1 as outcome. Constraint of k=6 is applied to all learners."

sl.perf <- y1_vaccine %>%
  filter(Learner == "SL") %>%
  mutate(AUCstr = ifelse(AUC %in% tail(sort(AUC), 1), paste0(AUCstr, "*"), AUCstr)) %>% 
  select(file, AUCstr) %>%
  mutate(file = fct_relevel(file,
                            c("1_only_matvars", 
                              "2_varset_EIA.log10d0", "3_varset_EIA.log10d14", "4_varset_EIA.log10d14overd0", "5_varset_EIA.log10cord",
                              "6_varset_PCA.log10d0", "7_varset_PCA.log10d14", "8_varset_PCA.log10d14overd0", "9_varset_PCA.log10cord",
                              "10_varset_RSVA.log10d0", "11_varset_RSVA.log10d14", "12_varset_RSVA.log10d14overd0", "13_varset_RSVA.log10cord",
                              "14_varset_RSVB.log10d0", "15_varset_RSVB.log10d14", "16_varset_RSVB.log10d14overd0", "17_varset_RSVB.log10cord",
                              "18_varset_d0_all4",
                              "19_varset_d14_all4", 
                              "20_varset_d14overd0_all4", "21_varset_cord_all4",
                              "22_varset_EIA_all4", "23_varset_PCA_all4", "24_varset_RSVA_all4", "25_varset_RSVB_all4",
                              "26_varset_d0_d14", "27_varset_d0_foldrise", "28_varset_d0_cord", "29_varset_all16"))) %>%
  arrange(file) %>%
  rename(`SL CV-AUC (95% CI)` = AUCstr,
         `Variable set` = file)


top.algo <- y1_vaccine %>%
  filter(AUCstr != "none") %>%
  group_by(file) %>%
  filter(AUC == max(AUC)) %>%
  ungroup() %>%
  mutate(AUCstr = ifelse(AUC %in% tail(sort(AUC), 1), paste0(AUCstr, "*"), AUCstr)) %>% 
  select(file, Learner, Screen, AUCstr) %>%
  rename(`Variable set` = file,
         `CV-AUC (95% CI)` = AUCstr) 

top.algo_none <- y1_vaccine %>%
  filter(AUCstr == "none") %>%
  distinct(file, .keep_all = TRUE) %>%
  ungroup() %>%
  select(file, Learner, Screen, AUCstr) %>%
  rename(`Variable set` = file,
         `CV-AUC (95% CI)` = AUCstr) 

top.algo <- bind_rows(top.algo, top.algo_none) %>%
  mutate(`Variable set` = fct_relevel(`Variable set`,
                                      c("1_only_matvars", 
                                        "2_varset_EIA.log10d0", "3_varset_EIA.log10d14", "4_varset_EIA.log10d14overd0", "5_varset_EIA.log10cord",
                                        "6_varset_PCA.log10d0", "7_varset_PCA.log10d14", "8_varset_PCA.log10d14overd0", "9_varset_PCA.log10cord",
                                        "10_varset_RSVA.log10d0", "11_varset_RSVA.log10d14", "12_varset_RSVA.log10d14overd0", "13_varset_RSVA.log10cord",
                                        "14_varset_RSVB.log10d0", "15_varset_RSVB.log10d14", "16_varset_RSVB.log10d14overd0", "17_varset_RSVB.log10cord",
                                        "18_varset_d0_all4","19_varset_d14_all4", "20_varset_d14overd0_all4", "21_varset_cord_all4",
                                        "22_varset_EIA_all4", "23_varset_PCA_all4", "24_varset_RSVA_all4", "25_varset_RSVB_all4",
                                        "26_varset_d0_d14", "27_varset_d0_foldrise", "28_varset_d0_cord", "29_varset_all16"))) %>%
  arrange(`Variable set`)

tab <- top.algo %>% full_join(sl.perf, by = "Variable set") %>%
  select(`Variable set`, `SL CV-AUC (95% CI)`, everything()) 

tab %>% write.csv(paste0(here(), "/correlates_report/objective3/input/SLperformance-vacc-y1.csv"))

## ----SLperformance-plac-y1, warning=kable_warnings--------------------------------------------------------------------------------------------------------------------
caption <- "Performance of Superlearner and the top-performing learner-screen combinations (CV-AUCs with 95\\% CIs) for each of the 29 variable sets using placebo group and endpoint 1 as outcome. Constraint of k=6 is applied to all learners."

sl.perf <- y1_placebo %>%
  filter(Learner == "SL") %>%
  mutate(AUCstr = ifelse(AUC %in% tail(sort(AUC), 1), paste0(AUCstr, "*"), AUCstr)) %>% 
  select(file, AUCstr) %>%
  mutate(file = fct_relevel(file, 
                            c("1_only_matvars",  
                              "2_varset_EIA.log10d0", "3_varset_EIA.log10d14", "4_varset_EIA.log10d14overd0", "5_varset_EIA.log10cord",
                              "6_varset_PCA.log10d0", "7_varset_PCA.log10d14", "8_varset_PCA.log10d14overd0", "9_varset_PCA.log10cord", 
                              "10_varset_RSVA.log10d0", "11_varset_RSVA.log10d14", "12_varset_RSVA.log10d14overd0", "13_varset_RSVA.log10cord",
                              "14_varset_RSVB.log10d0", "15_varset_RSVB.log10d14", "16_varset_RSVB.log10d14overd0", "17_varset_RSVB.log10cord", 
                              "18_varset_d0_all4","19_varset_d14_all4", "20_varset_d14overd0_all4", "21_varset_cord_all4",
                              "22_varset_EIA_all4", "23_varset_PCA_all4", "24_varset_RSVA_all4", "25_varset_RSVB_all4",
                              "26_varset_d0_d14", "27_varset_d0_foldrise", "28_varset_d0_cord", "29_varset_all16"))) %>%
  arrange(file) %>%
  rename(`SL CV-AUC (95% CI)` = AUCstr, 
         `Variable set` = file)


top.algo <- y1_placebo %>%
  group_by(file) %>%
  filter(AUC == max(AUC)) %>%
  ungroup() %>%
  mutate(AUCstr = ifelse(AUC %in% tail(sort(AUC), 1), paste0(AUCstr, "*"), AUCstr)) %>% 
  select(file, Learner, Screen, AUCstr) %>%
  rename(`Variable set` = file,
         `CV-AUC (95% CI)` = AUCstr) %>%
  mutate(`Variable set` = fct_relevel(`Variable set`,
                                      c("1_only_matvars", 
                                        "2_varset_EIA.log10d0", "3_varset_EIA.log10d14", "4_varset_EIA.log10d14overd0", "5_varset_EIA.log10cord",
                                        "6_varset_PCA.log10d0", "7_varset_PCA.log10d14", "8_varset_PCA.log10d14overd0", "9_varset_PCA.log10cord",
                                        "10_varset_RSVA.log10d0", "11_varset_RSVA.log10d14", "12_varset_RSVA.log10d14overd0", "13_varset_RSVA.log10cord",
                                        "14_varset_RSVB.log10d0", "15_varset_RSVB.log10d14", "16_varset_RSVB.log10d14overd0", "17_varset_RSVB.log10cord",
                                        "18_varset_d0_all4","19_varset_d14_all4", "20_varset_d14overd0_all4", "21_varset_cord_all4",
                                        "22_varset_EIA_all4", "23_varset_PCA_all4", "24_varset_RSVA_all4", "25_varset_RSVB_all4",
                                        "26_varset_d0_d14", "27_varset_d0_foldrise", "28_varset_d0_cord",
                                        "29_varset_all16"))) %>%
  arrange(`Variable set`)

tab <- top.algo %>% full_join(sl.perf, by = "Variable set") %>%
  select(`Variable set`, `SL CV-AUC (95% CI)`, everything()) 

tab %>% write.csv(paste0(here(), "/correlates_report/objective3/input/SLperformance-plac-y1.csv"))

## ----SLperformance-vacc-y2, warning=kable_warnings--------------------------------------------------------------------------------------------------------------------
caption <- "Performance of Superlearner and the top-performing learner-screen combinations (CV-AUCs with 95\\% CIs) for each of the 29 variable sets using vaccine group and endpoint 2 as outcome. Constraint of k=6 is applied to all learners."

sl.perf <- y2_vaccine %>%
  filter(Learner == "SL") %>%
  mutate(AUCstr = ifelse(AUC %in% tail(sort(AUC), 1), paste0(AUCstr, "*"), AUCstr)) %>% 
  select(file, AUCstr) %>%
  mutate(file = fct_relevel(file,
                            c("1_only_matvars", 
                              "2_varset_EIA.log10d0", "3_varset_EIA.log10d14", "4_varset_EIA.log10d14overd0", "5_varset_EIA.log10cord",
                              "6_varset_PCA.log10d0", "7_varset_PCA.log10d14", "8_varset_PCA.log10d14overd0", "9_varset_PCA.log10cord",
                              "10_varset_RSVA.log10d0", "11_varset_RSVA.log10d14", "12_varset_RSVA.log10d14overd0", "13_varset_RSVA.log10cord",
                              "14_varset_RSVB.log10d0", "15_varset_RSVB.log10d14", "16_varset_RSVB.log10d14overd0", "17_varset_RSVB.log10cord",
                              "18_varset_d0_all4",
                              "19_varset_d14_all4", 
                              "20_varset_d14overd0_all4", "21_varset_cord_all4",
                              "22_varset_EIA_all4", "23_varset_PCA_all4", "24_varset_RSVA_all4", "25_varset_RSVB_all4",
                              "26_varset_d0_d14", "27_varset_d0_foldrise", "28_varset_d0_cord", "29_varset_all16"))) %>%
  arrange(file) %>%
  rename(`SL CV-AUC (95% CI)` = AUCstr,
         `Variable set` = file)


top.algo <- y2_vaccine %>%
  filter(AUCstr != "none") %>%
  group_by(file) %>%
  filter(AUC == max(AUC)) %>%
  ungroup() %>%
  mutate(AUCstr = ifelse(AUC %in% tail(sort(AUC), 1), paste0(AUCstr, "*"), AUCstr)) %>% 
  select(file, Learner, Screen, AUCstr) %>%
  rename(`Variable set` = file,
         `CV-AUC (95% CI)` = AUCstr) 

top.algo_none <- y2_vaccine %>%
  filter(AUCstr == "none") %>%
  distinct(file, .keep_all = TRUE) %>%
  ungroup() %>%
  select(file, Learner, Screen, AUCstr) %>%
  rename(`Variable set` = file,
         `CV-AUC (95% CI)` = AUCstr) 

top.algo <- bind_rows(top.algo, top.algo_none) %>%
  mutate(`Variable set` = fct_relevel(`Variable set`,
                                      c("1_only_matvars", 
                                        "2_varset_EIA.log10d0", "3_varset_EIA.log10d14", "4_varset_EIA.log10d14overd0", "5_varset_EIA.log10cord",
                                        "6_varset_PCA.log10d0", "7_varset_PCA.log10d14", "8_varset_PCA.log10d14overd0", "9_varset_PCA.log10cord",
                                        "10_varset_RSVA.log10d0", "11_varset_RSVA.log10d14", "12_varset_RSVA.log10d14overd0", "13_varset_RSVA.log10cord",
                                        "14_varset_RSVB.log10d0", "15_varset_RSVB.log10d14", "16_varset_RSVB.log10d14overd0", "17_varset_RSVB.log10cord",
                                        "18_varset_d0_all4","19_varset_d14_all4", "20_varset_d14overd0_all4", "21_varset_cord_all4",
                                        "22_varset_EIA_all4", "23_varset_PCA_all4", "24_varset_RSVA_all4", "25_varset_RSVB_all4",
                                        "26_varset_d0_d14", "27_varset_d0_foldrise", "28_varset_d0_cord", "29_varset_all16"))) %>%
  arrange(`Variable set`)

tab <- top.algo %>% full_join(sl.perf, by = "Variable set") %>%
  select(`Variable set`, `SL CV-AUC (95% CI)`, everything()) 

tab %>% write.csv(paste0(here(), "/correlates_report/objective3/input/SLperformance-vacc-y2.csv"))

## ----SLperformance-plac-y2, warning=kable_warnings--------------------------------------------------------------------------------------------------------------------
caption <- "Performance of Superlearner and the top-performing learner-screen combinations (CV-AUCs with 95\\% CIs) for each of the 29 variable sets using placebo group and endpoint 2 as outcome. Constraint of k=6 is applied to all learners."

sl.perf <- y2_placebo %>%
  filter(Learner == "SL") %>%
  mutate(AUCstr = ifelse(AUC %in% tail(sort(AUC), 1), paste0(AUCstr, "*"), AUCstr)) %>% 
  select(file, AUCstr) %>%
  mutate(file = fct_relevel(file, 
                            c("1_only_matvars",  
                              "2_varset_EIA.log10d0", "3_varset_EIA.log10d14", "4_varset_EIA.log10d14overd0", "5_varset_EIA.log10cord",
                              "6_varset_PCA.log10d0", "7_varset_PCA.log10d14", "8_varset_PCA.log10d14overd0", "9_varset_PCA.log10cord", 
                              "10_varset_RSVA.log10d0", "11_varset_RSVA.log10d14", "12_varset_RSVA.log10d14overd0", "13_varset_RSVA.log10cord",
                              "14_varset_RSVB.log10d0", "15_varset_RSVB.log10d14", "16_varset_RSVB.log10d14overd0", "17_varset_RSVB.log10cord", 
                              "18_varset_d0_all4","19_varset_d14_all4", "20_varset_d14overd0_all4", "21_varset_cord_all4",
                              "22_varset_EIA_all4", "23_varset_PCA_all4", "24_varset_RSVA_all4", "25_varset_RSVB_all4",
                              "26_varset_d0_d14", "27_varset_d0_foldrise", "28_varset_d0_cord", "29_varset_all16"))) %>%
  arrange(file) %>%
  rename(`SL CV-AUC (95% CI)` = AUCstr, 
         `Variable set` = file)


top.algo <- y2_placebo %>%
  group_by(file) %>%
  filter(AUC == max(AUC)) %>%
  ungroup() %>%
  mutate(AUCstr = ifelse(AUC %in% tail(sort(AUC), 1), paste0(AUCstr, "*"), AUCstr)) %>% 
  select(file, Learner, Screen, AUCstr) %>%
  rename(`Variable set` = file,
         `CV-AUC (95% CI)` = AUCstr) %>%
  mutate(`Variable set` = fct_relevel(`Variable set`,
                                      c("1_only_matvars", 
                                        "2_varset_EIA.log10d0", "3_varset_EIA.log10d14", "4_varset_EIA.log10d14overd0", "5_varset_EIA.log10cord",
                                        "6_varset_PCA.log10d0", "7_varset_PCA.log10d14", "8_varset_PCA.log10d14overd0", "9_varset_PCA.log10cord",
                                        "10_varset_RSVA.log10d0", "11_varset_RSVA.log10d14", "12_varset_RSVA.log10d14overd0", "13_varset_RSVA.log10cord",
                                        "14_varset_RSVB.log10d0", "15_varset_RSVB.log10d14", "16_varset_RSVB.log10d14overd0", "17_varset_RSVB.log10cord",
                                        "18_varset_d0_all4","19_varset_d14_all4", "20_varset_d14overd0_all4", "21_varset_cord_all4",
                                        "22_varset_EIA_all4", "23_varset_PCA_all4", "24_varset_RSVA_all4", "25_varset_RSVB_all4",
                                        "26_varset_d0_d14", "27_varset_d0_foldrise", "28_varset_d0_cord",
                                        "29_varset_all16"))) %>%
  arrange(`Variable set`)

tab <- top.algo %>% full_join(sl.perf, by = "Variable set") %>%
  select(`Variable set`, `SL CV-AUC (95% CI)`, everything()) 

tab %>% write.csv(paste0(here(), "/correlates_report/objective3/input/SLperformance-plac-y2.csv"))

##############################################################################################################################
##############################################################################################################################
# Forest plots for selected models
############## y1 endpoint ########################
# vaccine group
png(paste0(here(), "/correlates_report/objective3/input/forest_y1_vaccine_8th_varset.png"), width=1000, height=1100)
top_learner <- make_forest_plot(y1_vaccine %>% filter(file=="8_varset_PCA.log10d14overd0"))
grid.arrange(top_learner$top_learner_nms_plot, top_learner$top_learner_plot, ncol=2)

# plot_grid(top$top_learner_nms_plot, top$top_learner_plot,
#           align = "none",
#           ncol=2,
#           rel_widths = c(1/2, 1/2))
dev.off()

# placebo group
png(paste0(here(), "/correlates_report/objective3/input/forest_y1_placebo_4th_varset.png"), width=1000, height=1100)
top_learner <- make_forest_plot(y1_vaccine %>% filter(file=="4_varset_EIA.log10d14overd0"))
grid.arrange(top_learner$top_learner_nms_plot, top_learner$top_learner_plot, ncol=2)

# plot_grid(top$top_learner_nms_plot, top$top_learner_plot,
#           align = "none",
#           ncol=2,
#           rel_widths = c(1/2, 1/2))
dev.off()
############## y2 endpoint ########################
# vaccine group
png(paste0(here(), "/correlates_report/objective3/input/forest_y2_vaccine_4th_varset.png"), width=1000, height=1100)
top_learner <- make_forest_plot(y2_vaccine %>% filter(file=="4_varset_EIA.log10d14overd0"))
grid.arrange(top_learner$top_learner_nms_plot, top_learner$top_learner_plot, ncol=2)

# plot_grid(top$top_learner_nms_plot, top$top_learner_plot, 
#           align = "none", 
#           ncol=2,
#           rel_widths = c(1/2, 1/2))
dev.off()


# png(paste0(here(), "/correlates_report/objective3/input/forest_vaccine_22nd_varset.png"), width=1000, height=1100)
# top_learner <- make_forest_plot(y2_vaccine %>% filter(file=="22_varset_EIA_all4"))
# grid.arrange(top_learner$top_learner_nms_plot, top_learner$top_learner_plot, ncol=2)
# 
# # plot_grid(top$top_learner_nms_plot, top$top_learner_plot, 
# #           align = "none", 
# #           ncol=2,
# #           rel_widths = c(1/2, 1/2))
# dev.off()


# placebo group
png(paste0(here(), "/correlates_report/objective3/input/forest_y2_placebo_23rd_varset.png"), width=1000, height=1100)
top_learner <- make_forest_plot(y2_placebo %>% filter(file=="23_varset_PCA_all4"))
grid.arrange(top_learner$top_learner_nms_plot, top_learner$top_learner_plot, ncol=2)

# plot_grid(top$top_learner_nms_plot, top$top_learner_plot,
#           align = "none",
#           ncol=2,
#           rel_widths = c(1/2, 1/2))
dev.off()
# png(paste0(here(), "/correlates_report/objective3/input/forest_y2_placebo_8th_varset.png"), width=1000, height=1100)
# top_learner <- make_forest_plot(y2_placebo %>% filter(file=="8_varset_PCA.log10d14overd0"))
# grid.arrange(top_learner$top_learner_nms_plot, top_learner$top_learner_plot, ncol=2)
# 
# # plot_grid(top$top_learner_nms_plot, top$top_learner_plot, 
# #           align = "none", 
# #           ncol=2,
# #           rel_widths = c(1/2, 1/2))
# dev.off()



##############################################################################################################################
##############################################################################################################################
# Find predictors selected in a Superlearner run 
# y1 endpoint
# Vaccine group
load("N:/cavd/Objective 4/GH-VAP/ID33-Wilson/Data/adata/objective3/results_with_allvars_n6/y1_vaccine/SLFITS/fits_y1_vaccine_8_varset_PCA.log10d14overd0.rda")
fit$coef %>% 
  as.data.frame() %>%
  tibble::rownames_to_column(var = "Learners") %>%
  rename(`Weights` = ".") %>%
  write.csv(paste0(here(), "/correlates_report/objective3/input/y1_vaccine_8_varset_PCA.log10d14overd0_weights.csv"))

fit[["fitLibrary"]]$screen_univariate_logistic_pval_plus_exposure_SL.step_All$object$coefficients %>% 
  as.data.frame() %>%
  tibble::rownames_to_column(var = "Predictors") %>%
  rename(`Coefficient` = ".") %>%
  mutate(`Odds Ratio` = exp(`Coefficient`)) %>%
  write.csv(paste0(here(), "/correlates_report/objective3/input/y1_vaccine_8_varset_PCA.log10d14overd0_learner.csv"))
rm(fit)

# Placebo group
load("N:/cavd/Objective 4/GH-VAP/ID33-Wilson/Data/adata/objective3/results_with_allvars_n6/y1_placebo/SLFITS/fits_y1_placebo_4_varset_EIA.log10d14overd0.rda")
fit$coef %>% 
  as.data.frame() %>%
  tibble::rownames_to_column(var = "Learners") %>%
  rename(`Weights` = ".") %>%
  write.csv(paste0(here(), "/correlates_report/objective3/input/y1_placebo_4_varset_EIA.log10d14overd0_weights.csv"))

fit[["fitLibrary"]]$screen_univariate_logistic_pval_plus_exposure_SL.step_All$object$coefficients %>% 
  as.data.frame() %>%
  tibble::rownames_to_column(var = "Predictors") %>%
  rename(`Coefficient` = ".") %>%
  mutate(`Odds Ratio` = exp(`Coefficient`)) %>%
  write.csv(paste0(here(), "/correlates_report/objective3/input/y1_placebo_4_varset_EIA.log10d14overd0_learner.csv"))
rm(fit)


# y2 endpoint
# Vaccine group
load("N:/cavd/Objective 4/GH-VAP/ID33-Wilson/Data/adata/objective3/results_with_allvars_n6/y2_vaccine/SLFITS/fits_y2_vaccine_4_varset_EIA.log10d14overd0.rda")
fit[["fitLibrary"]]$screen_all_plus_exposure_SL.glm_All$object$coefficients %>% 
  as.data.frame() %>%
  tibble::rownames_to_column(var = "Predictors") %>%
  rename(`Coefficient` = ".") %>%
  mutate(`Odds Ratio` = exp(`Coefficient`)) %>%
  write.csv(paste0(here(), "/correlates_report/objective3/input/coef_y2_vaccine_4_varset_EIA.log10d14overd0.csv"))
rm(fit)
# load("N:/cavd/Objective 4/GH-VAP/ID33-Wilson/Data/adata/objective3/results_with_allvars_n6/y2_vaccine/SLFITS/fits_y2_vaccine_22_varset_EIA_all4.rda")
# fit[["fitLibrary"]]$screen_all_plus_exposure_SL.glm_All$object$coefficients %>% 
#   as.data.frame() %>%
#   tibble::rownames_to_column(var = "Predictors") %>%
#   rename(`Coefficient` = ".") %>%
#   mutate(`Odds Ratio` = exp(`Coefficient`)) %>%
#   write.csv("H:/RSVcorrelatesAnalysis/correlates_report/objective3/input/coef_y2_vaccine_22_varset_EIA_all4.csv")
# rm(fit)
# Placebo group
load("N:/cavd/Objective 4/GH-VAP/ID33-Wilson/Data/adata/objective3/results_with_allvars_n6/y2_placebo/SLFITS/fits_y2_placebo_23_varset_PCA_all4.rda")
fit[["fitLibrary"]]$screen_glmnet_plus_exposure_SL.bayesglm_All$object$coefficients %>% 
  as.data.frame() %>%
  tibble::rownames_to_column(var = "Predictors") %>%
  rename(`Coefficient` = ".") %>%
  mutate(`Odds Ratio` = exp(`Coefficient`)) %>%
  write.csv(paste0(here(), "/correlates_report/objective3/input/coef_y2_placebo_23_varset_PCA_all4.csv"))
rm(fit)
# load("N:/cavd/Objective 4/GH-VAP/ID33-Wilson/Data/adata/objective3/results_with_allvars_n6/y2_placebo/SLFITS/fits_y2_placebo_8_varset_PCA.log10d14overd0.rda")
# library(glmnet)
# coef(fit[["fitLibrary"]]$screen_glmnet_plus_exposure_SL.glmnet_All$object) %>% 
#   as.matrix() %>% 
#   as.data.frame() %>%
#   tibble::rownames_to_column(var = "Predictors") %>%
#   rename(`Coefficient` = "1") %>%
#   mutate(`Odds Ratio` = exp(`Coefficient`)) %>%
#   write.csv("H:/RSVcorrelatesAnalysis/correlates_report/objective3/input/coef_y2_placebo_8_varset_PCA.log10d14overd0.csv")
# rm(fit)

######################################################################################################################
######################################################################################################################
# Appendix section 
# Make forest plots 
varset = relevel_fileColumn(y1_vaccine) %>% .$file %>% unique()
# print(varset) # make sure order is correct

for (i in c("y1", "y2")) {
  for (j in 1:29) { 
    # print(i)
    # print(j)
    # print(match(j, varset))
    if(i == "y1"){
      top <- make_forest_plot_2panels(y1_vaccine %>% filter(file == varset[j]),
                                      y1_placebo %>% filter(file == varset[j]))
      png(paste0(here(), "/correlates_report/objective3/input/y1_forest_varset_", j, ".png"), width=1000, height=1100)
      print(plot_grid(top$top_learner_nms_plot_vacc, top$top_learner_plot_vacc, top$top_learner_nms_plot_plac, top$top_learner_plot_plac,
                align = "none",
                ncol=4,
                rel_widths = c(1/2, 1/2, 1/5, 1/2)))
      dev.off()
    }
    if(i == "y2"){
      top <- make_forest_plot_2panels(y2_vaccine %>% filter(file == varset[j]),
                                      y2_placebo %>% filter(file == varset[j]))
      png(paste0(here(), "/correlates_report/objective3/input/y2_forest_varset_", j, ".png"), width=1000, height=1100)
      print(plot_grid(top$top_learner_nms_plot_vacc, top$top_learner_plot_vacc, top$top_learner_nms_plot_plac, top$top_learner_plot_plac,
                align = "none",
                ncol=4,
                rel_widths = c(1/2, 1/2, 1/5, 1/2)))
      dev.off()
    }
  }
}

######################################################################################################################
  

