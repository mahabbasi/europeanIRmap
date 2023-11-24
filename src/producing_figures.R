# ----------------------------- Function to produce Figures--------------------------------------

# create the variable important plot for both steps ---> Figure 7. ------------
ggvimp <- function(in_rftuned_step1, in_rftuned_step2, in_predvars, varnum = 23, spatial_rsp=FALSE) {
  
  rsmp_res_step1 <- get_outerrsmp(in_rftuned_step1, spatial_rsp=spatial_rsp)
  rsmp_res_step2 <- get_outerrsmp(in_rftuned_step2, spatial_rsp=spatial_rsp)
  #Get variable importance and format them
  varimp_basic_step1 <- weighted_vimportance_nestedrf(rfresamp = rsmp_res_step1,
                                                      pvalue = FALSE) %>%
    merge(., in_predvars, by.x='varnames', by.y = "varname") %>%
    .[, `:=`(varnames = factor(varnames, varnames[order(-imp_wmean)]),
             Category = factor(Category,
                               levels = c('Anthropogenic', 'Climate', 'Geology', 'Hydrology', 
                                          'Lakes','Landcover', 'Physiography'))
    )] %>%
    setorder(-imp_wmean)
  
  varimp_basic_step2 <- weighted_vimportance_nestedrf(rfresamp = rsmp_res_step2,
                                                      pvalue = FALSE) %>%
    merge(., in_predvars, by.x='varnames', by.y = "varname") %>%
    .[, `:=`(varnames = factor(varnames, varnames[order(-imp_wmean)]),
             Category = factor(Category,
                               levels = c('Anthropogenic', 'Climate', 'Geology', 'Hydrology', 
                                          'Lakes','Landcover', 'Physiography'))
    )] %>%
    setorder(-imp_wmean)
  #Plot step 1 RF
  outp_step1 <- varimp_basic_step1[1:dim(in_predvars)[1]] %>% 
    mutate(abbrevation = fct_reorder(abbrevation, desc(imp_wmean))) %>% 
    ggplot(. ,aes(x=abbrevation,
                  color =Category, fill=Category)) +
    geom_bar(aes(y=imp_wmean), stat = 'identity', alpha=0.7) +
    geom_errorbar(aes(ymin=imp_wmean-imp_wsd, ymax=imp_wmean+imp_wsd)) +
    scale_x_discrete(labels = function(x) {
      stringr::str_wrap(x, width = 27)
    },
    limits=rev) +
    scale_fill_manual(values=c('#FFCCFF', '#fdb462', '#696868','#3399FF', '#80b1d3','#b3de69', '#bc80bd'),
                      drop=FALSE) +
    scale_color_manual(values=c('#FFCCFF', '#fdb462', '#696868','#3399FF', '#80b1d3','#b3de69', '#bc80bd'),
                       drop=FALSE) +
    theme_classic(18) +
    theme(axis.text = element_text(color = 'black')) +
    scale_y_continuous(expand=c(0,0), position = 'left') +
    labs(y = "", x='', title = '              Step 1 RF')+
    coord_flip(ylim=c(0, max(varimp_basic_step1[, max(imp_wmean+imp_wsd)+10], 100))) 
  #Plot step 2
  outp_step2 <- varimp_basic_step2[1:dim(in_predvars)[1]] %>% 
    mutate(abbrevation = fct_reorder(abbrevation, desc(imp_wmean))) %>% 
    ggplot(. ,aes(x=abbrevation,
                  color =Category, fill=Category)) +
    geom_bar(aes(y=imp_wmean),position = position_dodge(), stat = 'identity', alpha=0.7) +
    
    geom_errorbar(aes(ymin=imp_wmean-imp_wsd, ymax=imp_wmean+imp_wsd)) +
    scale_y_reverse () +
    # scale_x_discrete() +
    coord_flip () +
    scale_x_discrete(labels = function(x) {
      stringr::str_wrap(x, width = 27)
    },
    limits=rev, position = "top") +
    scale_fill_manual(values=c('#FFCCFF', '#fdb462', '#696868','#3399FF', '#80b1d3','#b3de69', '#bc80bd'),
                      drop=FALSE) +
    scale_color_manual(values=c('#FFCCFF', '#fdb462', '#696868','#3399FF', '#80b1d3','#b3de69', '#bc80bd'),
                       drop=FALSE) +
    theme_classic(18) +
    theme(legend.position = 'none', axis.text = element_text(color = 'black')) +
    labs(y = "", x = '', title = '      Step 2 RF')
  
  pout <- cowplot::plot_grid(outp_step1, outp_step2, rel_widths = c(2, 1), label_size = 18) +
    cowplot::draw_label("Variable importance (actual impurity reduction)", 
                        x=0.5, y=  0, vjust=-0.5, angle= 0, size = 18) +
    cowplot::draw_label("Predictor", x=  0, y=0.5, vjust= 1.5, angle=90, size = 18) +
    cowplot::draw_label("Predictor", x=  1, y=0.5, vjust= -1.5, angle=90, size = 18)
  
  return(pout)
}

# Create a data table of predictors' name and category
predname_df <- function(in_task, in_task_all = NULL, feat_name_vec = NULL, category = NULL){
  
  if (inherits(in_task, "Task")) {
    feat_name <- in_task$feature_names
    if (length(list(in_task_all)) == 1) {
      feat_name_all <- in_task_all$feature_names
    } else {
      feat_name_all <- in_task_all[[1]]$feature_names
    }
    
  } else
    feat_name <- feat_name_vec
  
  feat_name_abb <- c('Q', 'P_to_PET_ratio', 'Q_iav_cv', 'dor', 'drainage_area', 'glacier_frac',
                     'gwr_to_runoff_ratio', 'irri_frac_im', 'irri_frac', 'karst_frac', 'karst_status', 'land_cover',
                     'lake_frac', 'Q_mean_p12', 'Q_mean_p3', 'Q_min_p12', 'Q_min_p3', 'pot_nat_vegetation',
                     'pop_dens_im', 'pop_dens', 'Q_iav_sd', 'slope', 'wet_days')
  if (is.null(category)) {
    cat <- c("Hydrology", "Climate", "Hydrology", "Anthropogenic", "Physiography", "Landcover", "Hydrology",
             "Anthropogenic", "Anthropogenic", "Geology", "Geology",
             "Landcover", "Lakes",  "Hydrology", "Hydrology", "Hydrology", "Hydrology",
             "Landcover", "Anthropogenic", "Anthropogenic", "Hydrology", "Physiography", "Climate")
  } else
    cat <- category
  
  predvars <- cbind(feat_name, feat_name_abb, cat) %>%
    `colnames<-`(c("varname", 'abbrevation',"Category")) %>% 
    as.data.table()
  
  predvars_out <- predvars[varname %in% feat_name]
  
  return(predvars_out)
  
}
# extract_pd_nestedrf -----------
extract_pd_nestedrf <- function(learner_id=1, in_rftuned, datdf,
                                selcols, nvariate, ngrid=20) {
  in_mod <- as.data.table(in_rftuned)[eval(learner_id),] #Go through data.table format to have access to both tasks and learners
  
  #Get fold-specific performance measure
  foldperf <- extract_impperf_nestedrf(in_rflearner = in_mod$learner,
                                       in_task = in_mod$task,
                                       imp=F, perf=T, pvalue=F)
  
  # selcols <- in_vimp_plot$data %>% #Can use that if extracting from tunredrf is expensive
  #   setorder(-imp_wmean) %>%
  #   .[colnums, variable]
  
  
  if (inherits(in_mod$learner[[1]]$learner, "GraphLearner")) {
    in_fit <- in_mod$learner[[1]]$learner$model$classif.ranger$model
  } else {
    in_fit <- in_mod$learner[[1]]$learner$model
  }
  
  ngridvec <- c(ngrid, ngrid)
  
  #Make dataset of all combinations of selected column names, two at a time
  if (nvariate == 1) {
    pdcomb <- lapply(selcols, function(i) {
      print(i)
      pdout <- edarf::partial_dependence(fit = in_fit, vars = c(i),
                                         n = ngridvec, data = datdf) %>% #Warning: does not work with data_table
        setDT %>%
        .[,(names(foldperf)) := foldperf] %>%
        .[, `:=`(var1=i)] %>%
        setnames(i, 'value1')
    }
    ) %>%
      do.call(rbind, .)
    
  } else if (nvariate == 2) {
    vargrid <- combn(selcols, 2, simplify=F) %>%
      do.call(rbind, .)
    
    #Get marginal distribution of the effect of two columns at a time
    pdcomb <- mapply(function(i, j) {
      pdout <- edarf::partial_dependence(fit = in_fit, vars = c(i, j),
                                         n = ngridvec,
                                         interaction = TRUE, data = datdf) %>% #Warning: does not work with data_table
        setDT %>%
        .[,(names(foldperf)) := foldperf] %>%
        .[, `:=`(var1=i, var2=j)] %>%
        setnames(c(i,j), c('value1', 'value2'))
      
      
      return(pdout)
    }, vargrid[,1], vargrid[,2], SIMPLIFY = FALSE) %>%
      do.call(rbind, .)
  } else {
    print('Warning: function cannot yet work with more than two variables at a time')
  }
  
  return(pdcomb)
}

# ggpartialdep ---------- 
ggpartialdep <- function (in_rftuned, in_predvars, colnums, ngrid, nodupli=T,
                          nvariate = 2, model_data,
                          parallel=T, spatial_rsp=FALSE) {
  
  in_predvars$unites <- c('(m3 sec-1 km-2)', '(-)', '(-)', '(% * 10)', '(km2)',
                          '(%)', '(-)', '(% * 100)', '(% * 100)', '(%)', '(-)', '(-)', '(% * 100)',
                          '(m3 sec-1 km-2)',
                          '(m3 sec-1 km-2)', '(m3 sec-1 km-2)', '(m3 sec-1 km-2)', '(-)', '(People per km2)',
                          '(People per km2)', '(m3 sec-1 km-2)', '(deg/100)', '(days/month * 10000)')
  
  in_predvars[, combined := paste(abbrevation, unites, sep = " ")]
  model_data[, gwr_to_runoff_ratio := gwr_to_runoff_ratio/100][
    ,wet_days := wet_days/100
    ][
      ,dor_pc_pva := dor_pc_pva/1000
      ][, lka_pc_use := lka_pc_use/100][, ire_pc_cse := ire_pc_cse/100][,ire_pc_use:=ire_pc_use/100]
  
  #Get outer resampling of interest
  rsmp_res <- get_outerrsmp(in_rftuned, spatial_rsp=spatial_rsp)
  
  #Get partial dependence across all folds and repeats
  nlearners <-with(rsmp_res$resampling$param_set$values, folds*repeats)
  datdf <- as.data.frame(model_data) #This may be shortened
  varimp <- weighted_vimportance_nestedrf(rsmp_res, pvalue=FALSE) %>%
    setorder(-imp_wmean)
  
  if (length(colnums) > nrow(varimp)) {
    colnums <- colnums[1:nrow(varimp)]
    cat('colnums argument exceeded the number of variables,
        reduced it to ', nrow(varimp), ' variables')
  }
  
  if (nodupli) {
    selcols <- as.character(
      varimp$varnames[!duplicated(substr(varimp$varnames, 1,3))][colnums])
  } else {
    selcols <- as.character(
      varimp$varnames[colnums])
  }
  
  if (parallel) {
    print(paste("Computing partial dependence with future.apply across", nlearners,
                "CV folds"))
    pd <- future.apply::future_lapply(seq_len(nlearners),
                                      extract_pd_nestedrf,
                                      in_rftuned = rsmp_res,
                                      datdf = datdf,
                                      selcols = selcols,
                                      nvariate = nvariate,
                                      ngrid = ngrid,
                                      future.scheduling = structure(TRUE,ordering = "random"),
                                      future.packages = c("data.table","edarf","ranger"))
    
  } else {
    print(paste("Computing partial dependence iteratively across", nlearners,
                "CV folds"))
    pd <- lapply(seq_len(nlearners),
                 extract_pd_nestedrf,
                 in_rftuned = rsmp_res,
                 datdf = datdf,
                 selcols = selcols,
                 nvariate = nvariate,
                 ngrid = ngrid)
  }
  
  #Get weighted mean
  varvec <- paste0('var', 1:nvariate)
  valvec <- paste0('value', 1:nvariate)
  
  pdformat <- do.call(rbind, pd) %>%
    setDT %>%
    .[, list(mean1 = weighted.mean(`1`, classif.bacc)),
      by= c(varvec, valvec)] %>%
    .[, variables := var1] %>%
    merge(., in_predvars, by.x='var1', by.y='varname')
  
  datdf2 <- as.data.table(datdf)[, target_class := as.numeric(as.character(target))]
  not_log_variables <- c('ai','gwr_to_runoff_ratio', 'pot_nat_vegetation', 'land_cover', 'cv', 'glacier_fraction',
                         'karst_fraction', 'karst_status', 'lka_pc_use', 'ppd_pk_cav', 'ppd_pk_uav', 'ire_pc_cse',
                         'ire_pc_use', 'slope', 'wet_days')
  if (nvariate ==1) {
    
    tileplots_list1 <- pdformat[!(var1 %in% not_log_variables),
                                list(list(ggplotGrob(
                                  ggplot(.SD, aes(x=value1, y=mean1)) +
                                    geom_line() +
                                    geom_rug(data=datdf2,
                                             aes_string(x=eval(var1),y='target_class'),
                                             alpha=1/3) +
                                    scale_y_continuous(name='Partial dependence (probability of intermittency)',
                                                       limits= c(min(mean1)-0.01, max(mean1)+0.01),  #c(0.25, 0.425),
                                                       expand=c(0,0))+
                                    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                                                  labels = trans_format("log10", math_format(10^.x)),
                                                  name=stringr::str_wrap(eval(combined), width = 30)) +
                                    theme_classic() +
                                    theme(text = element_text(size=16),
                                          axis.title.y = element_blank())
                                ))), by=.(var1)]
    
    tileplots_list2 <- pdformat[(var1 %in% not_log_variables),
                                list(list(ggplotGrob(
                                  ggplot(.SD, aes(x=value1, y=mean1)) +
                                    geom_line() +
                                    geom_rug(data=datdf2,
                                             aes_string(x=eval(var1),y='target_class'),
                                             alpha=1/3) +
                                    scale_y_continuous(name='Partial dependence (probability of intermittency)',
                                                       limits= c(min(mean1)-0.01, max(mean1)+0.01),  #c(0.25, 0.425),
                                                       expand=c(0,0)) +
                                    scale_x_continuous(name=stringr::str_wrap(eval(combined), width = 30)) +
                                    theme_classic() +
                                    theme(text = element_text(size=16),
                                          axis.title.y = element_blank())
                                ))), by=.(var1)]
    
    tileplots_l <- rbind(tileplots_list1, tileplots_list2)
    setkey(tileplots_l, var1)
    tileplots_l <- tileplots_l[varimp$varnames]
    
  }
  
  pagelayout <-   lapply(1:(nrow(tileplots_l) %/% 9), function(p_i) {
    (p_i-1)*9+(1:9)
  })
  if (nrow(tileplots_l) %% 9 > 0) {
    pagelayout[[nrow(tileplots_l) %/% 9 + 1]] <-
      (nrow(tileplots_l) %/% 9)*9+(1:(nrow(tileplots_l) %% 9))
  }
  
  
  tileplots_multipl <- lapply(pagelayout, function(page) {
    print(page)
    return(do.call("grid.arrange",list(
      grobs=(tileplots_l[page,V1]),
      left = textGrob('Probability that the station-month is intermittent',
                      gp = gpar(fontface = "bold", fontsize = 16), rot = 90))))
  })
  return(tileplots_multipl)
}


### compute logNSE and NSE metrics of downscaled waterGAP
compute_logNSE_NSE_dswatergap <- function(observed_streamflow, ds_watergap, shp){
  
  observed_streamflow_daily <- observed_streamflow
  
  
  observed_streamflow_daily[, year_month := format(dates, "%Y-%m")]
  
  target_dt <- observed_streamflow_daily[, .(obs_month = mean(value)), by = .(gaugeid, year_month)][
    , dates := as.Date(paste0(year_month, "-01"))][
      ,.(gaugeid, dates, obs_month)
      ]
  
  dswaterGAP_streamflow <- ds_watergap[,.(gaugeid, dates, Q)]
  
  merged_dt <- dplyr::left_join(target_dt, dswaterGAP_streamflow, by = c('gaugeid', 'dates')) %>% 
    setnames(c('obs_month', 'Q'), c('obs', 'sim')) %>% 
    .[obs == 0, obs := 0.0000001] %>% .[sim == 0, sim := 0.0000001]
  
  gauging_ids = merged_dt[,gaugeid] %>% unique()
  lognse_dt <- lapply(seq_along(gauging_ids), function(i){
    selected_gaugeid = gauging_ids[i]
    selected_dt <- merged_dt[gaugeid == selected_gaugeid]
    bamr::logNSE(selected_dt$sim, selected_dt$obs)
  }) %>% 
    do.call('rbind',.) %>% 
    as.data.table() %>%
    .[,gaugeid := gauging_ids] %>% 
    setnames('V1', "logNSE")
  
  nse_dt <- lapply(seq_along(gauging_ids), function(i){
    selected_gaugeid = gauging_ids[i]
    selected_dt <- merged_dt[gaugeid == selected_gaugeid]
    bamr::NSE(selected_dt$sim, selected_dt$obs)
  }) %>% 
    do.call('rbind',.) %>% 
    as.data.table() %>%
    .[,gaugeid := gauging_ids] %>% 
    setnames('V1', "NSE")
  
  final_gauges <- shp %>% 
    dplyr::rename(gaugeid = 'gauge_d')
  
  first_join <- merge(final_gauges, nse_dt, by = 'gaugeid')
  output_func <-  merge(first_join, lognse_dt, by = 'gaugeid')
  
  return(output_func)
}

### produce the boxplot of logNSE and NSE performance metrics --> Figure 2
ggboxplot_nse_lognse <- function(in_data){
  
  #Set the label for the plot grouped by upstream area of gauging stations
  data_label_lognse = c("(0-2] \n(17/3)", "(2-5] \n(22/2)",
                        "(5-10] \n(52/4)", "(10-50] \n(395/65)",
                        "(50-500] \n(1786/242)",
                        "(500-2,500] \n(790/143)", "(2500-10,000] \n(366/91)", "> 10,000 \n(281/67)")
  data_label_nse = c("(0-2] \n(17/3)", "(2-5] \n(22/3)",
                     "(5-10] \n(52/4)", "(10-50] \n(395/48)",
                     "(50-500] \n(1786/262)",
                     "(500-2,500] \n(790/144)", "(2500-10,000] \n(366/82)", "> 10,000 \n(281/49)")
  
  #Categorize the gauging stations based on their upstream area
  in_data <- in_data %>% 
    sf::st_drop_geometry() %>% 
    setDT() %>% 
    .[,.(up_mdfd, logNSE, NSE)] %>%  .[,bins_lognse:=cut(up_mdfd, breaks = c(-Inf, 2, 5, 10, 50, 
                                                                             500, 2500, 10000, Inf),
                                                         labels = data_label_lognse)] %>% 
    .[,bins_nse:=cut(up_mdfd, breaks = c(-Inf, 2, 5, 10, 50, 
                                         500, 2500, 10000, Inf),
                     labels = data_label_nse)] 
  
  
  # Function to calculate boxplot percentiles
  bp.pctiles = function (x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95)) {
    r <- quantile(x, probs = probs, na.rm = TRUE)
    names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
    r
  }
  
  #Plot the boxplot of logNSE performance metric in 8 different categories
  outp_lognse <- in_data %>%
    .[,.(logNSE, bins_lognse)] %>%
    reshape2::melt(., id.vars = c("bins_lognse")) %>%
    ggplot(aes(bins_lognse, xend=bins_lognse, y = value))+
    geom_violin(color= 'blue', width = 0.5) +
    stat_summary(fun.data=bp.pctiles, geom="boxplot", width= 0.2) +
    scale_y_continuous(breaks=c(-1,-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1))+
    coord_cartesian(ylim = c(-1, 1)) +
    theme_bw() +
    labs(x = "Upstream area of streamflow gauging stations [Km2]
         (number of stations/number of stations not shown)",
         y = "NSE of Log streamflow")+
    scale_fill_discrete(name = expression(paste('Upstream area ', km^{2})),
                        labels = c('[0-2)', "(2-5]", "(5-10]", "(10-50]",
                                   "(50-500]", "(500-2,500]", "(2500-10,000]", "> 10,000"))+
    theme(legend.title = element_text(colour = "black",
                                      family = "Times New Roman",
                                      size =  14,
                                      face = "bold"),
          legend.text =  element_text(colour = "black",
                                      family = "Times New Roman",
                                      size =  11,
                                      face = "bold"),
          legend.position = "bottom",
          axis.title = element_text(colour = "black",
                                    family = "Times New Roman",
                                    size =  14,
                                    face = "bold"),
          axis.text = element_text(colour = "black",
                                   family = "Times New Roman",
                                   size =  13,
                                   face = "bold"),
          panel.grid.major = element_line(colour = "black", linewidth = 0.4),
          panel.grid.minor = element_line(colour = "black", linewidth = 0.2),
          axis.line = element_line(colour = 'black', linewidth = 2)
    )
  
  #Plot NSE boxplot
  outp_nse <- in_data %>%
    .[,.(NSE, bins_nse)] %>%
    reshape2::melt(., id.vars = c("bins_nse")) %>%
    ggplot(aes(bins_nse, xend=bins_nse, y = value))+
    geom_violin(color= 'blue', width = 0.5) + 
    stat_summary(fun.data=bp.pctiles, geom="boxplot", width= 0.2)+
    scale_y_continuous(breaks=c(-1,-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1))+
    coord_cartesian(ylim = c(-1, 1)) +
    theme_bw() +
    labs(x = "Upstream area of streamflow gauging stations [Km2]
         (number of stations/number of stations not shown)",
         y = "NSE")+
    scale_fill_discrete(name = expression(paste('Upstream area ', km^{2})),
                        labels = c('[0-2)', "(2-5]", "(5-10]", "(10-50]",
                                   "(50-500]", "(500-2,500]", "(2500-10,000]", "> 10,000"))+
    theme(legend.title = element_text(colour = "black",
                                      family = "Times New Roman",
                                      size =  14, 
                                      face = "bold"),
          legend.text =  element_text(colour = "black",
                                      family = "Times New Roman",
                                      size =  11, 
                                      face = "bold"),
          legend.position = "bottom",
          axis.title = element_text(colour = "black",
                                    family = "Times New Roman",
                                    size =  14,
                                    face = "bold"),
          axis.text = element_text(colour = "black",
                                   family = "Times New Roman",
                                   size =  13, 
                                   face = "bold"),
          panel.grid.major = element_line(colour = "black", linewidth = 0.4),
          panel.grid.minor = element_line(colour = "black", linewidth = 0.2),
          axis.line = element_line(colour = 'black', linewidth = 2)
    ) 
  
  ####### outputs ---------
  outp <- ggpubr::ggarrange(outp_nse, outp_lognse, ncol=2, common.legend = TRUE, legend="bottom")
  
  return(outp)
}

ggboxplot_intermittent_figS5 <- function(in_nse_data, in_model_data){
  
  intermittent_gaugeids <- in_model_data[target >0,] %>% .[,unique(gaugeid)]
  perennial_gaugeids <- setdiff(in_model_data[,unique(gaugeid)],intermittent_gaugeids)
  intermittent_data <- in_nse_data %>% filter(gaugeid %in% intermittent_gaugeids)
  #Set the label for the plot grouped by upstream area of gauging stations
  data_label_lognse = c("(0-2] \n(14/3)", "(2-5] \n(16/1)",
                        "(5-10] \n(33/0)", "(10-50] \n(143/21)",
                        "(50-500] \n(465/58)",
                        "(500-2,500] \n(136/27)", "(2500-10,000] \n(38/13)", "> 10,000 \n(42/26)")
  data_label_nse = c("(0-2] \n(14/2)", "(2-5] \n(16/3)",
                     "(5-10] \n(33/1)", "(10-50] \n(143/26)",
                     "(50-500] \n(465/113)",
                     "(500-2,500] \n(136/56)", "(2500-10,000] \n(38/18)", "> 10,000 \n(42/25)")
  
  #Categorize the gauging stations based on their upstream area
  intermittent_data <- intermittent_data %>% 
    sf::st_drop_geometry() %>% 
    setDT() %>% 
    .[,.(up_mdfd, logNSE, NSE)] %>%  .[,bins_lognse:=cut(up_mdfd, breaks = c(-Inf, 2, 5, 10, 50, 
                                                                             500, 2500, 10000, Inf),
                                                         labels = data_label_lognse)] %>% 
    .[,bins_nse:=cut(up_mdfd, breaks = c(-Inf, 2, 5, 10, 50, 
                                         500, 2500, 10000, Inf),
                     labels = data_label_nse)] 
  
  # Function to calculate boxplot percentiles
  bp.pctiles = function (x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95)) {
    r <- quantile(x, probs = probs, na.rm = TRUE)
    names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
    r
  }
  
  #Plot the boxplot of logNSE performance metric in 8 different categories
  outp_lognse <- intermittent_data %>%
    .[,.(logNSE, bins_lognse)] %>%
    reshape2::melt(., id.vars = c("bins_lognse")) %>%
    ggplot(aes(bins_lognse, xend=bins_lognse, y = value))+
    geom_violin(color= 'blue', width = 0.5) + 
    stat_summary(fun.data=bp.pctiles, geom="boxplot", width= 0.2)+
    scale_y_continuous(breaks=c(-1,-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1))+
    coord_cartesian(ylim = c(-1, 1)) +
    theme_bw() +
    labs(x = "Upstream area of streamflow gauging stations [Km2]
         (number of stations/number of stations not shown)",
         y = "NSE of Log streamflow")+
    scale_fill_discrete(name = expression(paste('Upstream area ', km^{2})),
                        labels = c('[0-2)', "(2-5]", "(5-10]", "(10-50]",
                                   "(50-500]", "(500-2,500]", "(2500-10,000]", "> 10,000"))+
    theme(legend.title = element_text(colour = "black",
                                      family = "Times New Roman",
                                      size =  14, 
                                      face = "bold"),
          legend.text =  element_text(colour = "black",
                                      family = "Times New Roman",
                                      size =  11, 
                                      face = "bold"),
          legend.position = "bottom",
          axis.title = element_text(colour = "black",
                                    family = "Times New Roman",
                                    size =  14,
                                    face = "bold"),
          axis.text = element_text(colour = "black",
                                   family = "Times New Roman",
                                   size =  13, 
                                   face = "bold"),
          panel.grid.major = element_line(colour = "black", linewidth = 0.4),
          panel.grid.minor = element_line(colour = "black", linewidth = 0.2),
          axis.line = element_line(colour = 'black', linewidth = 2)
    )
  #Plot NSE as a function of upstream area
  outp_nse <- intermittent_data %>%
    .[,.(NSE, bins_nse)] %>%
    reshape2::melt(., id.vars = c("bins_nse")) %>%
    ggplot(aes(bins_nse, xend=bins_nse, y = value))+
    geom_violin(color= 'blue', width = 0.5) + 
    stat_summary(fun.data=bp.pctiles, geom="boxplot", width= 0.2)+
    scale_y_continuous(breaks=c(-1,-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1))+
    coord_cartesian(ylim = c(-1, 1)) +
    theme_bw() +
    labs(x = "Upstream area of streamflow gauging stations [Km2]
         (number of stations/number of stations not shown)",
         y = "NSE")+
    scale_fill_discrete(name = expression(paste('Upstream area ', km^{2})),
                        labels = c('[0-2)', "(2-5]", "(5-10]", "(10-50]",
                                   "(50-500]", "(500-2,500]", "(2500-10,000]", "> 10,000"))+
    theme(legend.title = element_text(colour = "black",
                                      family = "Times New Roman",
                                      size =  14, 
                                      face = "bold"),
          legend.text =  element_text(colour = "black",
                                      family = "Times New Roman",
                                      size =  11, 
                                      face = "bold"),
          legend.position = "bottom",
          axis.title = element_text(colour = "black",
                                    family = "Times New Roman",
                                    size =  14,
                                    face = "bold"),
          axis.text = element_text(colour = "black",
                                   family = "Times New Roman",
                                   size =  13, 
                                   face = "bold"),
          panel.grid.major = element_line(colour = "black", linewidth = 0.4),
          panel.grid.minor = element_line(colour = "black", linewidth = 0.2),
          axis.line = element_line(colour = 'black', linewidth = 2)
    ) 
  
  outp <- ggpubr::ggarrange(outp_nse, outp_lognse, ncol=2, common.legend = TRUE, legend="bottom")
  return(outp)
}

ggboxplot_perennial_figS6 <- function(in_nse_data, in_model_data){
  
  intermittent_gaugeids <- in_model_data[target >0,] %>% .[,unique(gaugeid)]
  perennial_data <- in_nse_data %>% filter(!gaugeid %in% intermittent_gaugeids)
  #Set the label for the plot grouped by upstream area of gauging stations
  data_label_lognse = c("(0-2] \n(3/0)", "(2-5] \n(6/1)",
                        "(5-10] \n(19/4)", "(10-50] \n(252/44)",
                        "(50-500] \n(1321/184)",
                        "(500-2,500] \n(654/116)", "(2500-10,000] \n(328/78)", "> 10,000 \n(239/41)")
  data_label_nse = c("(0-2] \n(3/1)", "(2-5] \n(6/0)",
                     "(5-10] \n(19/3)", "(10-50] \n(252/22)",
                     "(50-500] \n(1321/149)",
                     "(500-2,500] \n(654/88)", "(2500-10,000] \n(328/64)", "> 10,000 \n(239/24)")
  
  #Categorize the gauging stations based on their upstream area
  perennial_data <- perennial_data %>% 
    sf::st_drop_geometry() %>% 
    setDT() %>% 
    .[,.(up_mdfd, logNSE, NSE)] %>%  .[,bins_lognse:=cut(up_mdfd, breaks = c(-Inf, 2, 5, 10, 50, 
                                                                             500, 2500, 10000, Inf),
                                                         labels = data_label_lognse)] %>% 
    .[,bins_nse:=cut(up_mdfd, breaks = c(-Inf, 2, 5, 10, 50, 
                                         500, 2500, 10000, Inf),
                     labels = data_label_nse)] 
  
  
  # Function to calculate boxplot percentiles
  bp.pctiles = function (x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95)) {
    r <- quantile(x, probs = probs, na.rm = TRUE)
    names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
    r
  }
  
  #Plot the boxplot of logNSE performance metric in 8 different categories
  outp_lognse <- perennial_data %>%
    .[,.(logNSE, bins_lognse)] %>%
    reshape2::melt(., id.vars = c("bins_lognse")) %>%
    ggplot(aes(bins_lognse, xend=bins_lognse, y = value))+
    geom_violin(color= 'blue', width = 0.5) + 
    stat_summary(fun.data=bp.pctiles, geom="boxplot", width= 0.2)+
    scale_y_continuous(breaks=c(-1,-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1))+
    coord_cartesian(ylim = c(-1, 1)) +
    theme_bw() +
    labs(x = "Upstream area of streamflow gauging stations [Km2]
         (number of stations/number of stations not shown)",
         y = "NSE of Log streamflow")+
    scale_fill_discrete(name = expression(paste('Upstream area ', km^{2})),
                        labels = c('[0-2)', "(2-5]", "(5-10]", "(10-50]",
                                   "(50-500]", "(500-2,500]", "(2500-10,000]", "> 10,000"))+
    theme(legend.title = element_text(colour = "black",
                                      family = "Times New Roman",
                                      size =  14, 
                                      face = "bold"),
          legend.text =  element_text(colour = "black",
                                      family = "Times New Roman",
                                      size =  11, 
                                      face = "bold"),
          legend.position = "bottom",
          axis.title = element_text(colour = "black",
                                    family = "Times New Roman",
                                    size =  14,
                                    face = "bold"),
          axis.text = element_text(colour = "black",
                                   family = "Times New Roman",
                                   size =  13, 
                                   face = "bold"),
          panel.grid.major = element_line(colour = "black", linewidth = 0.4),
          panel.grid.minor = element_line(colour = "black", linewidth = 0.2),
          axis.line = element_line(colour = 'black', linewidth = 2)
    )
  
  outp_nse <- perennial_data %>%
    .[,.(NSE, bins_nse)] %>%
    reshape2::melt(., id.vars = c("bins_nse")) %>%
    ggplot(aes(bins_nse, xend=bins_nse, y = value))+
    geom_violin(color= 'blue', width = 0.5) + 
    stat_summary(fun.data=bp.pctiles, geom="boxplot", width= 0.2)+
    scale_y_continuous(breaks=c(-1,-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1))+
    coord_cartesian(ylim = c(-1, 1)) +
    theme_bw() +
    labs(x = "Upstream area of streamflow gauging stations [Km2]
         (number of stations/number of stations not shown)",
         y = "NSE")+
    scale_fill_discrete(name = expression(paste('Upstream area ', km^{2})),
                        labels = c('[0-2)', "(2-5]", "(5-10]", "(10-50]",
                                   "(50-500]", "(500-2,500]", "(2500-10,000]", "> 10,000"))+
    theme(legend.title = element_text(colour = "black",
                                      family = "Times New Roman",
                                      size =  14, 
                                      face = "bold"),
          legend.text =  element_text(colour = "black",
                                      family = "Times New Roman",
                                      size =  11, 
                                      face = "bold"),
          legend.position = "bottom",
          axis.title = element_text(colour = "black",
                                    family = "Times New Roman",
                                    size =  14,
                                    face = "bold"),
          axis.text = element_text(colour = "black",
                                   family = "Times New Roman",
                                   size =  13, 
                                   face = "bold"),
          panel.grid.major = element_line(colour = "black", linewidth = 0.4),
          panel.grid.minor = element_line(colour = "black", linewidth = 0.2),
          axis.line = element_line(colour = 'black', linewidth = 2)
    ) 
  
  outp <- ggpubr::ggarrange(outp_nse, outp_lognse, ncol=2, common.legend = TRUE, legend="bottom")
  return(outp)
}

### compute the ratio of no-flow observed and predicted -> figure 4 a and b ----
compute_noflow_ratio <- function(in_model, gauges_shp, in_model_data, outdir){
  
  dt <- in_model$rf_outer$prediction() %>% as.data.table()
  dt <- dt %>% 
    .[,c('mean_prob0', 'mean_prob1') := .(mean(prob.0), mean(prob.1)),
      by = row_ids] %>% 
    distinct(., row_ids,.keep_all = TRUE) %>% .[order(row_ids)]
  
  dt[, response := ifelse(mean_prob0 >= 0.5, 0, 1)]
  
  dt[,c('gaugeid', 'dates', 'month_date') := .(in_model_data$gaugeid, in_model_data$dates,
                                               format(in_model_data$dates, '%Y-%m'))]
  
  dt[, c('count', 'no_flow_obs', 'no_flow_pre') := 
       .(.N, sum(truth == 1), sum(response == 1)), by = gaugeid]
  dt[, c('ratio_no_flow_obs', 'ratio_no_flow_pre') := 
       .((no_flow_obs / count) * 100, (no_flow_pre / count) * 100), by = .(gaugeid)][
         ,c('ratio') := .(ratio_no_flow_pre / ratio_no_flow_obs),  by = .(gaugeid)
         ]
  
  distincted_dt <- unique(dt, by = 'gaugeid') %>%
    .[, .(gaugeid, count, no_flow_obs, no_flow_pre, 
          ratio_no_flow_obs, ratio_no_flow_pre, ratio)]
  
  distincted_dt[is.na(ratio), ratio := 9999]
  
  final_gauges <- gauges_shp %>%
    rename(gaugeid = gauge_d)
  
  output <- left_join(final_gauges, distincted_dt, by = 'gaugeid')
  output %>% sf::st_write(., dsn=file.path(here::here(), outdir, 'gaugingstations_figure4ab_step1.shp'))
  return(output)
}
## compute the correlation and NSE for model step1 ->> figure 4c -----
compute_corr_nse <- function(in_model, gauges_shp, in_model_data, outdir){
  
  dt <- in_model$rf_outer$prediction() %>% as.data.table()
  
  dt <- dt %>% .[,c('mean_prob0', 'mean_prob1') := .(mean(prob.0), mean(prob.1)), by = row_ids] %>% 
    distinct(., row_ids,.keep_all = TRUE) %>% .[order(row_ids)]
  dt[, response := ifelse(mean_prob0 >= 0.5, 0, 1)]
  dt[,c('gaugeid', 'dates', 'date_year') := .(in_model_data$gaugeid, in_model_data$dates,
                                              format(in_model_data$dates, '%Y'))]
  
  dt[, response := factor(response)]
  dt[, c('num_obs', 'num_pre') := .(sum(as.numeric(levels(truth)[truth])),
                                    sum(as.numeric(levels(response)[response]))),
     by = .(date_year, gaugeid)]
  
  dt[, c('corr', 'nse') := .(cor(num_pre, num_obs), NSE(num_pre, num_obs)), by = gaugeid]
  distincted_dt <- unique(dt, by = 'gaugeid') %>%
    .[, .(gaugeid, corr, nse)]
  distincted_dt[is.na(corr), corr := 9999]
  distincted_dt[is.na(nse), nse := 9999]
  final_gauges <- gauges_shp %>% rename(gaugeid = gauge_d)
  output <- left_join(final_gauges, distincted_dt, by = 'gaugeid')
  output %>% sf::st_write(., dsn=file.path(here::here(), outdir, 'corr_nse_step1.shp'))
  return(output)
}

### Making plots final version ------
#### plot boxplot figure 3.a and b. 
ggboxplot_lognse_fig3 <- function(in_data){
  
  #Set the label for the plot grouped by upstream area of gauging stations
  data_label_lognse = c("(0-2] \n(17/3)", "(2-5] \n(22/2)",
                        "(5-10] \n(52/4)", "(10-50] \n(395/65)",
                        "(50-500] \n(1786/242)",
                        "(500-2,500] \n(790/143)", "(2500-10,000] \n(366/91)", "> 10,000 \n(281/67)")
  data_label_nse = c("(0-2] \n(17/3)", "(2-5] \n(22/3)",
                     "(5-10] \n(52/4)", "(10-50] \n(395/48)",
                     "(50-500] \n(1786/262)",
                     "(500-2,500] \n(790/144)", "(2500-10,000] \n(366/82)", "> 10,000 \n(281/49)")
  
  #Categorize the gauging stations based on their upstream area
  in_data <- in_data %>% 
    sf::st_drop_geometry() %>% 
    setDT() %>% 
    .[,.(up_mdfd, logNSE, NSE)] %>%  .[,bins_lognse:=cut(up_mdfd, breaks = c(-Inf, 2, 5, 10, 50, 
                                                                             500, 2500, 10000, Inf),
                                                         labels = data_label_lognse)] %>% 
    .[,bins_nse:=cut(up_mdfd, breaks = c(-Inf, 2, 5, 10, 50, 
                                         500, 2500, 10000, Inf),
                     labels = data_label_nse)] 
  
  
  # Function to calculate boxplot percentiles
  bp.pctiles = function (x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95)) {
    r <- quantile(x, probs = probs, na.rm = TRUE)
    names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
    r
  }
  
  #Plot the boxplot of logNSE performance metric in 8 different categories
  
  outp_lognse <- in_data %>%
    .[,.(logNSE, bins_lognse)] %>%
    reshape2::melt(., id.vars = c("bins_lognse")) %>%
    ggplot(aes(bins_lognse, xend=bins_lognse, y = value))+
    geom_violin(color= 'blue', width = 0.5) +
    stat_summary(fun.data=bp.pctiles, geom="boxplot", width= 0.2, linewidth=0.8) +
    scale_y_continuous(breaks=c(-1,-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1))+
    coord_flip(ylim = c(-1, 1)) +
    theme_bw(18) +
    labs(x = "",
         y = "NSE of Log streamflow")+
    scale_fill_discrete(name = expression(paste('Upstream area ', km^{2})),
                        labels = c('[0-2)', "(2-5]", "(5-10]", "(10-50]",
                                   "(50-500]", "(500-2,500]", "(2500-10,000]", "> 10,000"))+
    theme(axis.text = element_text(colour = 'black', size = 16),
          legend.position = "bottom",
          panel.grid.major = element_line(colour = "black", linewidth = 0.3),
          panel.grid.minor = element_line(colour = "black", linewidth = 0.1),
          axis.line = element_line(colour = 'black', linewidth = 2)
    )
  
  outp_nse <- in_data %>%
    .[,.(NSE, bins_nse)] %>%
    reshape2::melt(., id.vars = c("bins_nse")) %>%
    ggplot(aes(bins_nse, xend=bins_nse, y = value))+
    geom_violin(color= 'blue', width = 0.5) + 
    stat_summary(fun.data=bp.pctiles, geom="boxplot", width= 0.2, linewidth=0.8)+
    scale_y_continuous(breaks=c(-1,-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1)) +
    coord_flip(ylim = c(-1, 1)) +
    # coord_cartesian(ylim = c(-1, 1)) 
    theme_bw(18) +
    labs(x = "Upstream area of streamflow gauging stations [Km2]
         (number of stations/number of stations not shown)",
         y = "NSE")+
    scale_fill_discrete(name = expression(paste('Upstream area ', km^{2})),
                        labels = c('[0-2)', "(2-5]", "(5-10]", "(10-50]",
                                   "(50-500]", "(500-2,500]", "(2500-10,000]", "> 10,000"))+
    theme(axis.text = element_text(colour = 'black', size = 16),
          legend.position = "bottom",
          panel.grid.major = element_line(colour = "black", linewidth = 0.3),
          panel.grid.minor = element_line(colour = "black", linewidth = 0.1),
          axis.line = element_line(colour = 'black', linewidth = 2)
    ) 
  
  ####### outputs ---------
  outp <- ggpubr::ggarrange(outp_nse, outp_lognse, ncol=2, common.legend = TRUE, legend="bottom")

  return(outp)
}

ggboxplot_intermittent_figS5 <- function(in_nse_data, in_model_data){
  
  intermittent_gaugeids <- in_model_data[target >0,] %>% .[,unique(gaugeid)]
  perennial_gaugeids <- setdiff(in_model_data[,unique(gaugeid)],intermittent_gaugeids)
  intermittent_data <- in_nse_data %>% filter(gaugeid %in% intermittent_gaugeids)
  #Set the label for the plot grouped by upstream area of gauging stations
  data_label_lognse = c("(0-2] \n(14/3)", "(2-5] \n(16/1)",
                        "(5-10] \n(33/0)", "(10-50] \n(143/21)",
                        "(50-500] \n(465/58)",
                        "(500-2,500] \n(136/27)", "(2500-10,000] \n(38/13)", "> 10,000 \n(42/26)")
  data_label_nse = c("(0-2] \n(14/2)", "(2-5] \n(16/3)",
                     "(5-10] \n(33/1)", "(10-50] \n(143/26)",
                     "(50-500] \n(465/113)",
                     "(500-2,500] \n(136/56)", "(2500-10,000] \n(38/18)", "> 10,000 \n(42/25)")
  
  #Categorize the gauging stations based on their upstream area
  intermittent_data <- intermittent_data %>% 
    sf::st_drop_geometry() %>% 
    setDT() %>% 
    .[,.(up_mdfd, logNSE, NSE)] %>%  .[,bins_lognse:=cut(up_mdfd, breaks = c(-Inf, 2, 5, 10, 50, 
                                                                             500, 2500, 10000, Inf),
                                                         labels = data_label_lognse)] %>% 
    .[,bins_nse:=cut(up_mdfd, breaks = c(-Inf, 2, 5, 10, 50, 
                                         500, 2500, 10000, Inf),
                     labels = data_label_nse)] 
  
  # Function to calculate boxplot percentiles
  bp.pctiles = function (x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95)) {
    r <- quantile(x, probs = probs, na.rm = TRUE)
    names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
    r
  }
  
  
  #Plot the boxplot of logNSE performance metric in 8 different categories
  outp_lognse <- intermittent_data %>%
    .[,.(logNSE, bins_lognse)] %>%
    reshape2::melt(., id.vars = c("bins_lognse")) %>%
    ggplot(aes(bins_lognse, xend=bins_lognse, y = value))+
    geom_violin(color= 'blue', width = 0.5) +
    stat_summary(fun.data=bp.pctiles, geom="boxplot", width= 0.2, linewidth=0.8) +
    scale_y_continuous(breaks=c(-1,-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1))+
    coord_flip(ylim = c(-1, 1)) +
    theme_bw(18) +
    labs(x = "",
         y = "NSE of Log streamflow")+
    scale_fill_discrete(name = expression(paste('Upstream area ', km^{2})),
                        labels = c('[0-2)', "(2-5]", "(5-10]", "(10-50]",
                                   "(50-500]", "(500-2,500]", "(2500-10,000]", "> 10,000"))+
    theme(axis.text = element_text(colour = 'black', size = 16),
          legend.position = "bottom",
          panel.grid.major = element_line(colour = "black", linewidth = 0.3),
          panel.grid.minor = element_line(colour = "black", linewidth = 0.1),
          axis.line = element_line(colour = 'black', linewidth = 2)
    )
  
  outp_nse <- intermittent_data %>%
    .[,.(NSE, bins_nse)] %>%
    reshape2::melt(., id.vars = c("bins_nse")) %>%
    ggplot(aes(bins_nse, xend=bins_nse, y = value))+
    geom_violin(color= 'blue', width = 0.5) + 
    stat_summary(fun.data=bp.pctiles, geom="boxplot", width= 0.2, linewidth=0.8)+
    scale_y_continuous(breaks=c(-1,-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1)) +
    coord_flip(ylim = c(-1, 1)) +
    # coord_cartesian(ylim = c(-1, 1)) 
    theme_bw(18) +
    labs(x = "Upstream area of streamflow gauging stations [Km2]
         (number of stations/number of stations not shown)",
         y = "NSE")+
    scale_fill_discrete(name = expression(paste('Upstream area ', km^{2})),
                        labels = c('[0-2)', "(2-5]", "(5-10]", "(10-50]",
                                   "(50-500]", "(500-2,500]", "(2500-10,000]", "> 10,000"))+
    theme(axis.text = element_text(colour = 'black', size = 16),
          legend.position = "bottom",
          panel.grid.major = element_line(colour = "black", linewidth = 0.3),
          panel.grid.minor = element_line(colour = "black", linewidth = 0.1),
          axis.line = element_line(colour = 'black', linewidth = 2)
    ) 
  
  outp <- ggpubr::ggarrange(outp_nse, outp_lognse, ncol=2, common.legend = TRUE, legend="bottom")
  return(outp)
}


ggboxplot_perennial_figS6 <- function(in_nse_data, in_model_data){
  
  intermittent_gaugeids <- in_model_data[target >0,] %>% .[,unique(gaugeid)]
  perennial_data <- in_nse_data %>% filter(!gaugeid %in% intermittent_gaugeids)
  #Set the label for the plot grouped by upstream area of gauging stations
  data_label_lognse = c("(0-2] \n(3/0)", "(2-5] \n(6/1)",
                        "(5-10] \n(19/4)", "(10-50] \n(252/44)",
                        "(50-500] \n(1321/184)",
                        "(500-2,500] \n(654/116)", "(2500-10,000] \n(328/78)", "> 10,000 \n(239/41)")
  data_label_nse = c("(0-2] \n(3/1)", "(2-5] \n(6/0)",
                     "(5-10] \n(19/3)", "(10-50] \n(252/22)",
                     "(50-500] \n(1321/149)",
                     "(500-2,500] \n(654/88)", "(2500-10,000] \n(328/64)", "> 10,000 \n(239/24)")
  
  #Categorize the gauging stations based on their upstream area
  perennial_data <- perennial_data %>% 
    sf::st_drop_geometry() %>% 
    setDT() %>% 
    .[,.(up_mdfd, logNSE, NSE)] %>%  .[,bins_lognse:=cut(up_mdfd, breaks = c(-Inf, 2, 5, 10, 50, 
                                                                             500, 2500, 10000, Inf),
                                                         labels = data_label_lognse)] %>% 
    .[,bins_nse:=cut(up_mdfd, breaks = c(-Inf, 2, 5, 10, 50, 
                                         500, 2500, 10000, Inf),
                     labels = data_label_nse)] 
  
  
  # Function to calculate boxplot percentiles
  bp.pctiles = function (x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95)) {
    r <- quantile(x, probs = probs, na.rm = TRUE)
    names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
    r
  }
  
  #Plot the boxplot of logNSE performance metric in 8 different categories
  outp_lognse <- perennial_data %>%
    .[,.(logNSE, bins_lognse)] %>%
    reshape2::melt(., id.vars = c("bins_lognse")) %>%
    ggplot(aes(bins_lognse, xend=bins_lognse, y = value))+
    geom_violin(color= 'blue', width = 0.5) +
    stat_summary(fun.data=bp.pctiles, geom="boxplot", width= 0.2, linewidth=0.8) +
    scale_y_continuous(breaks=c(-1,-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1))+
    coord_flip(ylim = c(-1, 1)) +
    theme_bw(18) +
    labs(x = "",
         y = "NSE of Log streamflow")+
    scale_fill_discrete(name = expression(paste('Upstream area ', km^{2})),
                        labels = c('[0-2)', "(2-5]", "(5-10]", "(10-50]",
                                   "(50-500]", "(500-2,500]", "(2500-10,000]", "> 10,000"))+
    theme(axis.text = element_text(colour = 'black', size = 16),
          legend.position = "bottom",
          panel.grid.major = element_line(colour = "black", linewidth = 0.3),
          panel.grid.minor = element_line(colour = "black", linewidth = 0.1),
          axis.line = element_line(colour = 'black', linewidth = 2)
    )
  
  outp_nse <- perennial_data %>%
    .[,.(NSE, bins_nse)] %>%
    reshape2::melt(., id.vars = c("bins_nse")) %>%
    ggplot(aes(bins_nse, xend=bins_nse, y = value))+
    geom_violin(color= 'blue', width = 0.5) + 
    stat_summary(fun.data=bp.pctiles, geom="boxplot", width= 0.2, linewidth=0.8)+
    scale_y_continuous(breaks=c(-1,-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1)) +
    coord_flip(ylim = c(-1, 1)) +
    # coord_cartesian(ylim = c(-1, 1)) 
    theme_bw(18) +
    labs(x = "Upstream area of streamflow gauging stations [Km2]
         (number of stations/number of stations not shown)",
         y = "NSE")+
    scale_fill_discrete(name = expression(paste('Upstream area ', km^{2})),
                        labels = c('[0-2)', "(2-5]", "(5-10]", "(10-50]",
                                   "(50-500]", "(500-2,500]", "(2500-10,000]", "> 10,000"))+
    theme(axis.text = element_text(colour = 'black', size = 16),
          legend.position = "bottom",
          panel.grid.major = element_line(colour = "black", linewidth = 0.3),
          panel.grid.minor = element_line(colour = "black", linewidth = 0.1),
          axis.line = element_line(colour = 'black', linewidth = 2)
    ) 
  
  outp <- ggpubr::ggarrange(outp_nse, outp_lognse, ncol=2, common.legend = TRUE, legend="bottom")
  return(outp)
}


ggplotConfusionMatrix_fig5 <- function(rftuned_step2, model_data){
  
  
  data_step2 <- model_data[target > 0]
  dt <- rftuned_step2$rf_outer$prediction() %>% as.data.table() %>% 
    .[,c('mean_prob1', 'mean_prob2', 'mean_prob3', 'mean_prob4') := .(mean(prob.1), mean(prob.2),
                                                                      mean(prob.3), mean(prob.4)),
      by = row_ids] %>% 
    distinct(., row_ids,.keep_all = TRUE) %>% .[order(row_ids)] %>% 
    .[, c('prob.1', 'prob.2', 'prob.3', 'prob.4') := NULL]
  
  max_values <- pmax(dt$mean_prob1, dt$mean_prob2, dt$mean_prob3, dt$mean_prob4)
  dt[, class_max := ifelse(max_values == mean_prob1, 1,
                           ifelse(max_values == mean_prob2, 2,
                                  ifelse(max_values == mean_prob3, 3, 4)))]
  
  dt[,class_max := as.factor(class_max)]
  dt[,c('gaugeid', 'dates', 'month_date') := .(data_step2$gaugeid, data_step2$dates,
                                               format(data_step2$dates, '%Y-%m'))]
  
  m <- caret::confusionMatrix(dt$class_max, dt$truth)
  
  data_c <-  mutate(group_by(as.data.frame(m$table), Reference ), percentage = 
                      percent(Freq/sum(Freq)))
  
  p <- ggplot(data = data_c,
              aes(x = Reference, y = Prediction)) +
    geom_tile(aes(fill = Freq), colour = "white") +
    scale_fill_gradient(low = "white", high = "orange") +
    geom_text(aes(x = Reference, y = Prediction, label = percentage), size = 5) +
    geom_text(aes(x = Reference, y = Prediction, label = Freq), vjust = -2, size = 5) +
    scale_x_discrete(labels=c('1-5 no-flow days', '6-15 no-flow days', '16-29 no-flow days', '30-31 no-flow days'))+
    scale_y_discrete(labels=c('1-5 no-flow days', '6-15 no-flow days', '16-29 no-flow days', '30-31 no-flow days')) +
    theme_bw(16) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 30, vjust = 0.5, family = 'TT Arial', color = 'black'),
          axis.text.y = element_text(angle = 30, vjust = 0.5, family = 'TT Arial', color = 'black')) +
    labs(x = 'Observed', y = 'Predicted')
  
  return(p)
}

compute_shp_fig6 <- function(in_model_step1, in_model_step2, in_model_data, gauges_shp,outdir){
  
  data_step2 <- in_model_data[target > 0] #filter intermittent station-months
  
  dt <- in_model_step2$rf_outer$prediction() %>% as.data.table() %>% 
    .[,c('mean_prob1', 'mean_prob2', 'mean_prob3', 'mean_prob4') := .(mean(prob.1), mean(prob.2),
                                                                      mean(prob.3), mean(prob.4)),
      by = row_ids] %>% 
    distinct(., row_ids,.keep_all = TRUE) %>% .[order(row_ids)] %>% 
    .[, c('prob.1', 'prob.2', 'prob.3', 'prob.4') := NULL]
  
  max_values <- pmax(dt$mean_prob1, dt$mean_prob2, dt$mean_prob3, dt$mean_prob4)
  dt[, class_max := ifelse(max_values == mean_prob1, 1,
                           ifelse(max_values == mean_prob2, 2,
                                  ifelse(max_values == mean_prob3, 3, 4)))]
  
  dt[,class_max := as.factor(class_max)]
  dt[,c('gaugeid', 'dates', 'month_date') := .(data_step2$gaugeid, data_step2$dates,
                                               format(data_step2$dates, '%Y-%m'))]
  
  
  new_df <- dt %>% 
    group_by(gaugeid) %>% 
    dplyr::summarise(count = n(),
                     ratio = sum(truth == 1 & class_max == 1 |
                                   truth == 2 & class_max == 2 |
                                   truth == 3 & class_max == 3 |
                                   truth == 4 & class_max == 4) / count * 100)
  
  final_gauges <- gauges_shp %>% 
    rename(gaugeid = gauge_d)
  
  merge_st_fig6a <- merge(new_df, final_gauges, by = "gaugeid") %>% 
    sf::st_as_sf()
  merge_st_fig6a %>% sf::write_sf(., dsn=file.path(outdir, "correctly_classfied_figure6a.shp"))
  
  cmean_cobs <- dt %>% 
    mutate(
      truth = as.numeric(levels(truth)[truth]),
      class_max = as.numeric(levels(class_max)[class_max])) %>%
    group_by(gaugeid) %>% 
    dplyr::summarise(cm_co = mean(class_max) - mean(truth))
  
  merge_st_fig6b <- merge(cmean_cobs, final_gauges, by = "gaugeid") %>% 
    sf::st_as_sf()
  merge_st_fig6b %>% sf::write_sf(., dsn=file.path(outdir, "cmean_cobs_step2_figure6b.shp"))
  
  cor_both <- dt %>% 
    mutate(
      truth = as.numeric(levels(truth)[truth]),
      class_max = as.numeric(levels(class_max)[class_max])) %>%
    group_by(gaugeid) %>% 
    dplyr::summarise(cor_pr = cor(truth, class_max, method = 'pearson'),
                     cor_spr = cor(truth, class_max, method = 'spearman'))
  
  cor_both <- cor_both %>% as.data.table() %>%
    .[is.na(cor_pr), cor_pr := -9999] %>% 
    .[is.na(cor_spr), cor_spr := -9999]
  
  merge_st_fig6c <- merge(cor_both, final_gauges, by = "gaugeid") %>% 
    sf::st_as_sf()
  merge_st_fig6c %>% sf::write_sf(., dsn=file.path(outdir, "correlation_both_step2_figure6c.shp"))
  
  # figure 6d - include the perennial station-months along with intermittent for 885 intermittent gauging stations
  dt_step1 <- in_model_step1$rf_outer$prediction() %>% as.data.table()
  
  intermittent_gaugeid <- dt[, base::unique(gaugeid)]
  
  dt_step1 <- dt_step1 %>% 
    .[,c('mean_prob0', 'mean_prob1') := .(mean(prob.0), mean(prob.1)), by = row_ids] %>% 
    distinct(., row_ids,.keep_all = TRUE) %>% .[order(row_ids)]
  
  dt_step1[, response := ifelse(mean_prob0 >= 0.5, 0, 1)]
  
  dt_step1[,c('gaugeid', 'dates', 'date_year') := .(in_model_data$gaugeid, in_model_data$dates,
                                                    format(in_model_data$dates, '%Y'))]
  
  dt_step1_interm <- dt_step1[gaugeid %in% intermittent_gaugeid][,response := as.factor(response)]
  dt_step1_interm[dt, response := i.class_max, on=c('gaugeid', 'dates')]
  dt_step1_interm[dt, truth := i.truth, on=c('gaugeid', 'dates')]
  
  cor_withzero <- dt_step1_interm %>% 
    mutate(
      truth = as.numeric(levels(truth)[truth]),
      response = as.numeric(levels(response)[response])) %>%
    group_by(gaugeid) %>% 
    dplyr::summarise(cor_pr = cor(truth, response, method = 'pearson'),
                     cor_spr = cor(truth, response, method = 'spearman'))
  
  cor_withzero <- cor_withzero %>% as.data.table() %>%
    .[is.na(cor_pr), cor_pr := -9999] %>%
    .[is.na(cor_spr), cor_spr := -9999]
  
  merge_st_fig6d <- merge(cor_withzero, final_gauges, by = "gaugeid") %>% 
    sf::st_as_sf()
  merge_st_fig6d %>% sf::write_sf(., dsn=file.path(outdir, "correlation_both_step2_withzero_figure6d.shp"))
  
  out_ls <- list(merge_st_fig6a, merge_st_fig6b, merge_st_fig6c, merge_st_fig6d)
  
  return(out_ls)
}

select_inter_month <- function(in_reachdt, path_static, path_reach_shp, outdir, in_date='date_2019_8'){
  
  # reach_data <- fst::read_fst(path = in_reachdt)
  
  net_eu_id_dt <-  sf::read_sf(dsn = path_reach_shp)
  
  upa <- fst::read_fst(path = path_static,
                       columns = c('DRYVER_RIVID',"drainage_area")) %>% 
    as.data.table() %>% rename(upa = drainage_area,
                               DRYVER_RIV = DRYVER_RIVID)
  
  years <- 1981:2019
  col_names <- vector("character")
  for(i in seq_along(years)){
    for(j in 1:12){
      col_names[(i-1)*12 + j] <- paste0("date_",years[i],"_",j)
    }
  }
  
  in_reachdt <- in_reachdt %>% 
    `colnames<-`(c("DRYVER_RIV", col_names))
  
  
  selected_mon_mr <- base::merge(in_reachdt[,c("DRYVER_RIV", in_date), with=FALSE], upa, by = "DRYVER_RIV")
  
  sel_mon_st <- merge(net_eu_id_dt, selected_mon_mr, by = "DRYVER_RIV") %>% sf::st_as_sf() 
  sf::write_sf(sel_mon_st, dsn=file.path(outdir, paste0("all_reaches_eu", in_date, '.shp')))
}

ggline_inter_ts <- function(in_reachdt, in_reach_shp=NULL, region=FALSE){
  
  if (region == TRUE) {
    region_sf <-  sf::read_sf(in_reach_shp)
    region_ids <- region_sf %>% 
      sf::st_drop_geometry() %>%
      pull(DRYVER_RIV)
    
    in_reachdt <- in_reachdt %>% 
      filter(V1 %in% region_ids)
    
    classes_reach_mon_df <- lapply(1:5, function(i){
      out <- colSums(in_reachdt == (i-1))/dim(in_reachdt)[1] *100
      out[-1]
    }) %>% 
      do.call("rbind", .)
  }else {
    classes_reach_mon_df <- lapply(1:5, function(i){
      out <- colSums(in_reachdt == (i-1))/dim(in_reachdt)[1] *100
      out[-1]
    }) %>% 
      do.call("rbind", .)
  }
  
  df1 <- classes_reach_mon_df %>% t() %>% as.data.table() %>% 
    cbind(time = seq.Date(as.Date("1981-01-01"),
                          as.Date("2019-12-31"), "month")) %>% 
    `colnames<-`(c("Perennial", "1-5 no-flow days",
                   "6-15 no-flow days", "16-29 no-flow days",
                   "30-31 no-flow days", "time")) %>% 
    data.table::melt(., id.vars = "time",
                     measure.vars = c("Perennial", "1-5 no-flow days",
                                      "6-15 no-flow days", "16-29 no-flow days",
                                      "30-31 no-flow days"))
  pout <- df1 %>%
    ggplot(., aes(x = time, y = value, color = variable)) +
    geom_line()+
    geom_point(size = 0.5, color='black')+
    theme_bw(16) + 
    facet_wrap("variable", scales = "free", ncol = 1) +
    scale_color_manual( values = c("deepskyblue4", "green",
                                   "yellow", "darkorange4", "red")) +
    labs(y = "Fraction of reaches over Europe (%)",
         x = "date (month)") +
    theme(panel.grid.major.y = element_line(colour = "black"),
          legend.position = "none",
          axis.title = element_text(colour = "black",
                                    size =  16,
                                    face = "bold"),
          axis.text = element_text(colour = "black",
                                   size =  14))
  
  return(pout)
}
