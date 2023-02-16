library(gridExtra)
library(grid)

choose_learners <- function(cvaucDAT, varSet){
  cvaucDAT %>% filter(file == varSet) %>% 
    filter(!Learner %in% c("SL", "Discrete SL")) %>%
    arrange(-AUC) %>% .[1:2,] %>%
    bind_rows(cvaucDAT %>% filter(file == varSet) %>% 
                filter(Learner %in% c("SL", "Discrete SL")))
}

get_predictions <- function(cv_fit, cvaucDAT, varSet){
  
  top3 <- choose_learners(cvaucDAT, varSet)
  
  predict <- cv_fit[["library.predict"]] %>% as.data.frame() %>%
    bind_cols(cv_fit[["discreteSL.predict"]] %>% as.data.frame() %>% `colnames<-`(c("Discrete SL"))) %>%
    bind_cols(cv_fit[["SL.predict"]] %>% as.data.frame() %>% `colnames<-`(c("Super Learner"))) %>%
    bind_cols(cv_fit[["Y"]] %>% as.data.frame() %>% `colnames<-`(c("Y"))) %>%
    gather("algo", "pred", -Y) %>%
    filter(algo %in% c(top3$Screen_fromRun, "Super Learner", "Discrete SL"))
  
  predict %>%
    left_join(top3 %>% select(Screen_fromRun, Learner, Screen, AUC), by = c("algo" = "Screen_fromRun")) %>%
    mutate(learnerScreen = paste0(Learner, "_", Screen),
           learnerScreen = ifelse(Learner %in% c("SL", "Discrete SL"), algo, learnerScreen),
           AUCchar = format(round(AUC, 3), nsmall=3),
           learnerScreen = paste0(learnerScreen, " (", AUCchar, ")"),
           learnerScreen =  reorder(learnerScreen, -AUC))
}



# Plot SL, discrete.SL and topRanking learner-screen combinations
plot_roc_curves <- function(predict, cvaucDAT, varSet) {
  
  top3 <- choose_learners(cvaucDAT, varSet)
  
  roc.obj <- predict %>%
    group_by(algo) %>%
    nest() %>%
    mutate(pred.obj = purrr::map(data, ~ ROCR::prediction(.x$pred, .x$Y)),
           perf.obj = purrr::map(pred.obj, ~ ROCR::performance(.x, "tpr", "fpr")),
           roc.dat = purrr::map(perf.obj, ~ tibble(xval = .x@x.values[[1]],
                                                   yval = .x@y.values[[1]]))) 
  
  roc.obj %>%
    unnest(roc.dat) %>%
    select(algo, xval, yval) %>%
    ungroup() %>%
    left_join(top3 %>% select(Screen_fromRun, Learner, Screen, AUC), by = c("algo" = "Screen_fromRun")) %>%
    mutate(learnerScreen = paste0(Learner, "_", Screen),
           learnerScreen = ifelse(Learner %in% c("SL", "Discrete SL"), algo, learnerScreen),
           AUCchar = format(round(AUC, 3), nsmall=3),
           learnerScreen = paste0(learnerScreen, " (", AUCchar, ")"),
           learnerScreen =  reorder(learnerScreen, -AUC)) %>%
    ggplot(aes(x=xval, y=yval, col=learnerScreen)) +
    geom_step(lwd=2) +
    theme(legend.position = "top", 
          legend.direction = "vertical",
          legend.box = "horizontal") +
    colorblindr::scale_color_OkabeIto(order = c(8, 1, 4, 3, 2)) +
    #scale_color_manual(values = cols) +
    labs(x = "Cross-Validated False Positive Rate", y = "Cross-Validated True Positive Rate", col = "Model (CV-AUC)") +
    geom_abline(intercept =0 , slope = 1)
}



plot_predicted_probabilities <- function(pred){
  pred %>% 
    mutate(Ychar = ifelse(Y==0, "Control", "Case")) %>%
    ggplot(aes(x=Ychar, y=pred, color=Ychar)) + 
    geom_boxplot() +
    geom_jitter(width = 0.2) +
    theme_bw() +
    #scale_color_manual(c("blue", "red"))
    colorblindr::scale_color_OkabeIto(order = c(5, 1)) +
    facet_grid(cols = vars(learnerScreen)) +
    labs(y = "CV estimated probability of RSV disease", x = "") +
    theme(legend.position = "none",
          strip.text.x = element_text(size = 11),
          axis.text = element_text(size = 12),
          axis.title.y = element_text(size = 14))
  # theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=-1)) +
}



# Read in the CV-AUC dataset to get the top 5 SL performers 
load(file = "N:/cavd/Objective 4/GH-VAP/ID33-Wilson/Data/adata/objective3/results_with_allvars_n6/y1_vaccine/y1_vaccine.rda")
load(file = "N:/cavd/Objective 4/GH-VAP/ID33-Wilson/Data/adata/objective3/results_with_allvars_n6/y1_placebo/y1_placebo.rda")
load(file = "N:/cavd/Objective 4/GH-VAP/ID33-Wilson/Data/adata/objective3/results_with_allvars_n6/y2_vaccine/y2_vaccine.rda")
load(file = "N:/cavd/Objective 4/GH-VAP/ID33-Wilson/Data/adata/objective3/results_with_allvars_n6/y2_placebo/y2_placebo.rda")
#################################################################################################################################
# For y1 vaccine and selected variable set and plot pred.Prob with SL, Discrete SL and top 2 best-performing individual Learners
y1_vaccine <- y1_vaccine %>% 
  mutate(Screen_fromRun = sapply(strsplit(Screen_fromRun,"_All"), `[`, 1),
         Screen_fromRun = paste0(Screen_fromRun, "_", Learner, "_All"),
         Screen_fromRun = ifelse(Learner == "SL", "Super Learner",
                                 ifelse(Learner == "Discrete SL", "Discrete SL", Screen_fromRun)))
top2_vacc <- y1_vaccine %>% 
  filter(Learner=="SL") %>%
  arrange(-AUC) %>%
  dplyr::slice(1:2)

load(file = paste0("N:/cavd/Objective 4/GH-VAP/ID33-Wilson/Data/adata/objective3/results_with_allvars_n6/y1_vaccine/CVFITS/slfits_y1_vaccine_", top2_vacc$file[1], ".rda"))
pred <- get_predictions(cvfits[[1]], cvaucDAT = y1_vaccine, varSet = top2_vacc$file[1])
png(file = paste0(here(), "/correlates_report/objective3/input/predProb_y1_vacc_", top2_vacc$file[1], ".png"), width=900, height=600)
#p1 <- plot_roc_curves(pred, cvaucDAT = y2_vaccine, varSet = top2_vacc$file[1]) 
p2 <- plot_predicted_probabilities(pred)
# vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
# grid.newpage()
# pushViewport(viewport(layout = grid.layout(8, 4)))
# print(p1, vp = vplayout(1:4, 2:3))
# print(p2, vp = vplayout(5:8, 1:4))
print(p2)
dev.off()
# rm(cvfits, pred, p1, p2)
rm(cvfits, pred, p2)
#################################################################################################################################
# For y1 placebo and selected variable set and plot pred.Prob with SL, Discrete SL and top 2 best-performing individual Learners
y1_placebo <- y1_placebo %>% 
  mutate(Screen_fromRun = sapply(strsplit(Screen_fromRun,"_All"), `[`, 1),
         Screen_fromRun = paste0(Screen_fromRun, "_", Learner, "_All"),
         Screen_fromRun = ifelse(Learner == "SL", "Super Learner",
                                 ifelse(Learner == "Discrete SL", "Discrete SL", Screen_fromRun)))
top2_plac <- y1_placebo %>% 
  filter(Learner=="SL") %>%
  arrange(-AUC) %>%
  dplyr::slice(1:2)

load(file = paste0("N:/cavd/Objective 4/GH-VAP/ID33-Wilson/Data/adata/objective3/results_with_allvars_n6/y1_placebo/CVFITS/slfits_y1_placebo_", top2_plac$file[1], ".rda"))
pred <- get_predictions(cvfits[[1]], cvaucDAT = y1_placebo, varSet = top2_plac$file[1])
png(file = paste0(here(), "/correlates_report/objective3/input/predProb_y1_plac_", top2_plac$file[1], ".png"), width=900, height=600)
p2 <- plot_predicted_probabilities(pred)
print(p2)
dev.off()
rm(cvfits, pred, p2)

#################################################################################################################################
# For y2 vaccine and selected variable set and plot pred.Prob with SL, Discrete SL and top 2 best-performing individual Learners
y2_vaccine <- y2_vaccine %>% 
  mutate(Screen_fromRun = sapply(strsplit(Screen_fromRun,"_All"), `[`, 1),
         Screen_fromRun = paste0(Screen_fromRun, "_", Learner, "_All"),
         Screen_fromRun = ifelse(Learner == "SL", "Super Learner",
                                 ifelse(Learner == "Discrete SL", "Discrete SL", Screen_fromRun)))

top2_vacc <- y2_vaccine %>% 
  arrange(-AUC) %>%
  dplyr::slice(1:2)

load(file = paste0("N:/cavd/Objective 4/GH-VAP/ID33-Wilson/Data/adata/objective3/results_with_allvars_n6/y2_vaccine/CVFITS/slfits_y2_vaccine_", top2_vacc$file[1], ".rda"))
pred <- get_predictions(cvfits[[1]], cvaucDAT = y2_vaccine, varSet = top2_vacc$file[1])
png(file = paste0(here(), "/correlates_report/objective3/input/predProb_y2_vacc_", top2_vacc$file[1], ".png"), width=900, height=600)
p2 <- plot_predicted_probabilities(pred)
print(p2)
dev.off()
rm(cvfits, pred, p2)

#################################################################################################################################
# For y2 placebo and selected variable set and plot pred.Prob with SL, Discrete SL and top 2 best-performing individual Learners
y2_placebo <- y2_placebo %>% 
  mutate(Screen_fromRun = sapply(strsplit(Screen_fromRun,"_All"), `[`, 1),
         Screen_fromRun = paste0(Screen_fromRun, "_", Learner, "_All"),
         Screen_fromRun = ifelse(Learner == "SL", "Super Learner",
                                 ifelse(Learner == "Discrete SL", "Discrete SL", Screen_fromRun)))

top2_plac <- y2_placebo %>% 
  arrange(-AUC) %>%
  dplyr::slice(1:2)

load(file = paste0("N:/cavd/Objective 4/GH-VAP/ID33-Wilson/Data/adata/objective3/results_with_allvars_n6/y2_placebo/CVFITS/slfits_y2_placebo_", top2_plac$file[1], ".rda"))
pred <- get_predictions(cvfits[[1]], cvaucDAT = y2_placebo, varSet = top2_plac$file[1])
png(file = paste0(here(), "/correlates_report/objective3/input/predProb_y2_plac_", top2_plac$file[1], ".png"), width=900, height=600)
p2 <- plot_predicted_probabilities(pred)
print(p2)
dev.off()
rm(cvfits, pred, p2)



