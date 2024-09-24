run_tmle3_outcome_rates <- function(Y, A, confounders){
  
  df <- data.frame(Y, A, confounders)
  node_list <- list(
    W = names(confounders),
    A = "A",
    Y = "Y"
  )
  
  
  processed <- process_missing(df, node_list)
  df <- processed$data
  node_list <- processed$node_list
  
  
  ## Define ATE's outcome levels
  ate_spec <- tmle_ATE(
    treatment_level = "1",
    control_level = "0"
  )
  
  
  ## Define learners
  lrn_mean <- Lrnr_mean$new()
  
  lrn_glmnet <- Lrnr_glmnet$new()
  lrn_earth <- Lrnr_earth$new()
  lrn_tree <- Lrnr_glmtree$new()
  lrn_tree2 <- Lrnr_rpart$new()
  lrn_gbm <- Lrnr_lightgbm$new()
  lrn_rf <- Lrnr_ranger$new()
  
  
  # define metalearners appropriate to data types
  ls_metalearner <- make_learner( # leaner for binary outcome
    Lrnr_nnls
  ) 
  mn_metalearner <- make_learner(  # learner for multinomial outcome
    Lrnr_solnp,
    metalearner_linear_multinomial,
    loss_loglik_multinomial
  )
  c_metalearner <- make_learner(  # learner for continuout outcome
    Lrnr_solnp,
    metalearner_linear,
    loss_squared_error
  )
  
  
  sl_Y <- Lrnr_sl$new(
    learners = list(lrn_mean, lrn_glmnet, lrn_tree2, lrn_rf),
    metalearner = c_metalearner
  )
  
  if(length(unique(A)) == 2){
    sl_A <- Lrnr_sl$new(
      learners = list(lrn_glmnet, lrn_tree2, lrn_rf),
      metalearner = ls_metalearner
    )
  } else {
    sl_A <- Lrnr_sl$new(
      learners = list(lrn_glmnet, lrn_tree2, lrn_rf),
      metalearner = mn_metalearner
    )
  }
  
  
  learner_list <- list(A = sl_A, Y = sl_Y)
  
  
  # Do it for each level in one step
  tmle_task <- ate_spec$make_tmle_task(df, node_list)
  
  initial_likelihood <- ate_spec$make_initial_likelihood(
    tmle_task,
    learner_list
  )
  
  tsm_spec <- tmle_TSM_all()
  
  targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood)
  all_tsm_params <- tsm_spec$make_params(tmle_task, targeted_likelihood)
  
  ate_params <- lapply(all_tsm_params[-1],
                       function(tsm_param){define_param(
                         Param_delta,
                         targeted_likelihood,
                         delta_param_RR,
                         list(all_tsm_params[[1]], tsm_param)
                       )})
  
  all_params <- c(all_tsm_params, ate_params)
  
  tmle_fit_multiparam_base_cfd <- fit_tmle3(
    tmle_task, targeted_likelihood, all_params,
    targeted_likelihood$updater
  )
  
  return(tmle_fit_multiparam_base_cfd)
}

make_box_plots <- function(df_box, y_vars, x){
  meltData <- melt(df_box %>% dplyr::select(y_vars, x), id.vars = x) 
  names(meltData)[1] = "x"
  if(x == 'age_cat'){
    meltData$x <- factor(meltData$x, levels = c('18-30', '31-45', '46-65', '>65')) 
  } else {
    meltData$x <- as.factor(meltData$x)
  }
  meltData$item = sub("^(n_|h_m_)", "", meltData$variable)
  pi <- ggplot(meltData, aes(x = x, y = value, color = x)) +
    geom_boxplot(position = position_dodge(width = 1)) +  
    facet_wrap(~ item, scales = 'free', nrow = 1) +
    theme_bw() +
    theme(legend.position = 'none',
          axis.text.x = element_text(angle = -75, vjust = 0.5, hjust=0),
          axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) +
    labs(color="x", x = "", y = "count")
  return(pi)
}


make_box_plots_no_outliers <- function(df_box, y_vars, x){
  meltData <- melt(df_box %>% dplyr::select(y_vars, x), id.vars = x) 
  names(meltData)[1] = "x"
  if(x == 'age_cat'){
    meltData$x <- factor(meltData$x, levels = c('18-30', '31-45', '46-65', '>65')) 
  } else {
    meltData$x <- as.factor(meltData$x)
  }
  meltData$item = sub("^(n_|h_m_)", "", meltData$variable)
  
  # Custom calculation for the whiskers
  whisker_data <- meltData %>%
    group_by(item, x) %>%
    summarize(lower = min(value),
              upper = max(value)) %>%
    ungroup()
  
  pi <- ggplot(meltData, aes(x = x, y = value, color = x)) +
    geom_boxplot(outlier.shape = NA, position = position_dodge(width = 1)) +
    geom_linerange(data = whisker_data, aes(x = x, ymin = lower, ymax = upper), color = "black") +
    facet_wrap(~ item, scales = 'free', nrow = 1) +
    theme_bw() +
    theme(legend.position = 'none',
          axis.text.x = element_text(angle = -75, vjust = 0.5, hjust=0),
          axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) +
    labs(color="x", x = "", y = "count")
  return(pi)
}



make_box_plots_tmle_rlts <- function(df_box, y_type, rotate_x_lab, x_lab, num_row = 1){
  pi <- ggplot(df_box %>% filter(Y_type == y_type), 
               aes(x = group)) +
    geom_errorbar(aes(ymin=lower, ymax=upper, color=group), width=0.7) +
    geom_point(aes(y=tmle_est, color = group))+
    facet_wrap(~ Y_name, scales = 'free', nrow = num_row) +
    labs(y = "TMLE estimated rates", x = x_lab) + # subtitle = 'TMLE estimation with other baseline variables and SOFA score as confounders'
    theme_bw() +
    theme(legend.position = 'none',
          plot.title = element_text(hjust = 0.5),  # plot.subtitle = element_text(hjust = 0.5)
          strip.text = element_text(size = 12),
          # strip.background = element_rect(fill = "#FFFFF0"),  # element_blank(),
          axis.text = element_text(size = 10.5))
  
  if (rotate_x_lab) {
    pi <- pi + theme(axis.text.x = element_text(angle = -75, vjust = 0.5, hjust=0))
  }
  
  return(pi)
}

make_box_plots_tmle_rlts_text <- function(df_box, y_type){
  
  df_tmp <- df_box %>% 
    filter(Y_type == y_type) %>%
    mutate(lower_2dgt = round(lower, 2),
           upper_2dgt = round(upper, 2),
           tmle_est_2dgt = round(tmle_est, 2))
  
  pi <- ggplot(df_tmp, aes(x = group)) +
    geom_errorbar(aes(ymin=lower, ymax=upper, color=group), width=0.7) +
    geom_point(aes(y=tmle_est, color = group)) +
    geom_text(aes(y=upper_2dgt, label=upper_2dgt), vjust=-0.5, color="black") +
    geom_text(aes(y=lower_2dgt, label=lower_2dgt), vjust=1.5, color="black") +
    geom_text(aes(y=tmle_est_2dgt, label=tmle_est_2dgt), vjust=1.5, color="black") +
    facet_wrap(~ Y_name, scales = 'free', nrow = 1) +
    labs(y="TMLE estimated rates") +
    theme_bw() +
    theme(legend.position = 'none',
          axis.text.x = element_text(angle = -75, vjust = 0.5, hjust=0),
          axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5))
  
  
  return(pi)
}


make_plots_tmle_results <- function(rlt_list, trt_name, orders,
                                    y_vars = NULL, y_type_1 = 'n_vital',
                                    fig_width = 17, fig_height = 3,
                                    num_row = 1) {
  
  # prep the data for plots
  tmle_fit_all <- do.call("rbind", rlt_list) %>% as.data.frame()
  
  tmle_fit_all_tsm <- tmle_fit_all %>%
    filter(type == 'TSM') %>%
    mutate(group = gsub(".*A=(.*)\\}.*", "\\1", param),
           Y_type = ifelse(Y %in% var_name_list$num_vitals, "n_vital", 
                           ifelse(Y %in% var_name_list$hours_missing_vitals, "h_m_vital",
                                  ifelse(Y %in% var_name_list$num_labs, "n_lab", 'other'))),
           Y_name = sub("^(n_|h_m_)", "", Y))
  
  if (trt_name == "Gender"){
    tmle_fit_all_tsm$group = ifelse(tmle_fit_all_tsm$group == "1", "Female", "Male")
  } else if (trt_name %in% c("Insuranced")) {
    tmle_fit_all_tsm$group = ifelse(tmle_fit_all_tsm$group == "1", "Yes", "No")
  } else if (trt_name == "Language") {
    tmle_fit_all_tsm$group[tmle_fit_all_tsm$group == "ENGL"] = "English"
  }
  tmle_fit_all_tsm$group = factor(tmle_fit_all_tsm$group,
                                  levels = orders)
  
  if ((!is.null(y_vars)) & all(y_vars != "all")) {
    tmle_fit_all_tsm <- tmle_fit_all_tsm |>
      dplyr::filter(Y_name %in% y_vars)
  }
  
  # making plots
  rotate_x_lab <- ifelse(max(nchar(orders)) > 7 | length(orders) > 6, TRUE, FALSE)
  
  if (is.null(y_vars)) {
    g_vs_hm <- make_box_plots_tmle_rlts(tmle_fit_all_tsm, 'h_m_vital', rotate_x_lab, trt_name) +
      labs(title = "Marginal Average Hours with Missing Vital Signs Measurements") 
    g_vs_hm
    
    g_vs_n <- make_box_plots_tmle_rlts(tmle_fit_all_tsm, 'n_vital', rotate_x_lab, trt_name) +
      labs(title = "Marginal Average Vital Sign Measurement Frequencies") 
    g_vs_n
    
    g_lab_n <- make_box_plots_tmle_rlts(tmle_fit_all_tsm, 'n_lab', rotate_x_lab, trt_name) +
      labs(title = "Marginal Average Laboratory Testing Frequencies") 
    g_lab_n
    
    ggsave(file=here(plotFolder, paste0(x, "_tmle_rlts_box_plot_vital_h_m.png")),
           g_vs_hm, width = fig_width, height = fig_height, dpi = 1200)
    
    ggsave(file=here(plotFolder, paste0(x, "_tmle_rlts_box_plot_vital_n.png")),
           g_vs_n, width = fig_width, height = fig_height, dpi = 1200)
    
    ggsave(file=here(plotFolder, paste0(x, "_tmle_rlts_box_plot_lab_n.png")),
           g_lab_n, width = fig_width, height = fig_height, dpi = 1200)
  } else {
    all_y_str <- paste(y_vars, collapse = "_")
    if (y_type_1 == "h_m_vital") {
      g_vs_hm <- make_box_plots_tmle_rlts(tmle_fit_all_tsm, 'h_m_vital',
                                          rotate_x_lab, trt_name, num_row) +
        labs(title = "Marginal Average Hours with Missing Vital Signs Measurements") 
      g_vs_hm
      ggsave(file=here(plotFolder, paste0(x, "_tmle_rlts_box_plot_vital_h_m_small_", all_y_str, ".png")),
             g_vs_hm, width = fig_width, height = fig_height, dpi = 1200)
    } else if (y_type_1 == "n_vital") {
      g_vs_n <- make_box_plots_tmle_rlts(tmle_fit_all_tsm, 'n_vital',
                                         rotate_x_lab, trt_name, num_row) +
        labs(title = "Marginal Average Vital Sign Measurement Frequencies") 
      g_vs_n
      ggsave(file=here(plotFolder, paste0(x, "_tmle_rlts_box_plot_vital_n_small_", all_y_str, ".png")),
             g_vs_n, width = fig_width, height = fig_height, dpi = 1200)
    } else {
      g_lab_n <- make_box_plots_tmle_rlts(tmle_fit_all_tsm, 'n_lab',
                                          rotate_x_lab, trt_name, num_row) +
        labs(title = "Marginal Average Laboratory Testing Frequencies") 
      g_lab_n
      ggsave(file=here(plotFolder, paste0(x, "_tmle_rlts_box_plot_lab_n_small_", all_y_str, ".png")),
             g_lab_n, width = fig_width, height = fig_height, dpi = 1200)
    }
  }
}
