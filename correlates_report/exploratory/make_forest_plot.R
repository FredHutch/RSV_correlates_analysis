
get_fancy_screen_names <- function(avgs){
  split_screen = avgs %>%
    mutate(prek = str_sub(Screen, end = -7),
           k = str_sub(Screen, -5, -5))

  whatk = unique(split_screen %>% filter(k!="") %>% .$k)

  returnDat <- split_screen %>%
    mutate(fancyScreen = case_when(prek == "screen_highcor_random_plus_exposure" ~ paste0("highcor_random_k=", k),
                                   prek == "screen_glmnet_plus_exposure" ~ paste0("glmnet_k=", k),
                                   prek == "screen_univariate_logistic_pval_plus_exposure" ~ paste0("univar_logistic_pval_k=", k),
                                   prek == "screen_all_plus_exposure" ~ paste0("all_k=", k),
                                   prek == "" ~ paste0("k=", whatk),
                                   TRUE ~ as.character(Screen)),
           k = whatk)

  return(returnDat)
}


make_forest_plot <- function(avgs){
  lowestXTick <- floor(min(avgs$ci_ll)*10)/10
  top_learner_plot <- ggplot() +
    geom_pointrange(avgs %>% mutate(LearnerScreen = fct_reorder(LearnerScreen, AUC, .desc = F)), mapping=aes(x=LearnerScreen, y=AUC, ymin=ci_ll, ymax=ci_ul), size = 0.35, color="blue", fill="blue", shape=20) +
    coord_flip() +
    scale_y_continuous(breaks = seq(lowestXTick, 1, 0.1), labels = seq(lowestXTick, 1, 0.1), limits = c(lowestXTick, 1)) +
    theme_bw() +
    labs(y = "CV-AUC [95% CI]", x = "") +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          plot.margin=unit(c(1,-0.15,1,-0.15),"cm"),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"))

  total_learnerScreen_combos = length(avgs$LearnerScreen)

  avgs_withCoord <- avgs %>%
    select(Learner, Screen, AUCstr) %>%
    gather("columnVal", "strDisplay") %>%
    mutate(xcoord = case_when(columnVal=="Learner" ~ 1,
                              columnVal=="Screen" ~ 1.5,
                              columnVal=="AUCstr" ~ 2),
           ycoord = rep(total_learnerScreen_combos:1, 3))

  top_learner_nms_plot <- ggplot(avgs_withCoord, aes(x = xcoord, y = ycoord, label = strDisplay)) +
    geom_text(hjust=1, vjust=0, size=3) +
    xlim(0.7,2) +
    theme(plot.margin=unit(c(0.8,-0.15,1.2,-0.15),"cm"),
          axis.line=element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 2, color = "white"),
          axis.ticks = element_blank(),
          axis.title = element_blank())

  return(list(top_learner_plot = top_learner_plot, top_learner_nms_plot = top_learner_nms_plot))
}
