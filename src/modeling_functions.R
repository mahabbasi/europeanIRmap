#### -------------------- Develope RF model Step ONE ----------------
## classify 
classify_target <- function(tbl,target_col, breaks, labels){
  target_class <- tbl %>%
    pull(target_col) %>%
    cut(breaks = breaks,
        labels = labels)
  return(target_class)
}

# creating task for different problem type ---
create_task_stepone <- function(in_tbl, problem_type = "classif", ...){
  colnames(in_tbl) <- make.names(names(in_tbl),
                                 unique = TRUE)
  
  in_tbl = in_tbl[target > 0, target:=1]
  
  if (!inherits(in_tbl[,target], "factor")) {
    in_tbl[, target := as.factor(target)]
  }
  
  if (problem_type == "regr"){
    task <- mlr3spatiotempcv::TaskRegrST$new(
      id = "regression_no_flow_days",
      backend = in_tbl,
      target = target_name,
      coordinate_names = c("X", "Y")
    )
  }else if (problem_type == "classif") {
    
    task <- mlr3spatiotempcv::TaskClassifST$new(
      id = "binary_class",
      backend = in_tbl,
      target = "target",
      coordinate_names = c("X", "Y"))
    
  } else {
    stop("the problem type have to be one of regr and classif types,
         please check the problem_type argument.\n")
  }
  
  return(task)
  }
#------ get_oversamp_ratio -----------------
#' Get oversample ratio
#'
#' Identify minority class and compute ratio between the number of observations
#' in the majority and the minority classes of a binary
#' \link[mlr3]{TaskClassif }.
#'
#' @param in_task binary \link[mlr3]{TaskClassif }
#'
#' @return named list with following items:
#' \describe{
#' \item{minoclass} - Value of minority class (e.g. '1' for intermittent rivers)
#' \item{ratio} - numeric ratio between number of items in majority and minority
#' class
#' }
#'
#' @examples
#' \dontrun{
#' in_dt <- data.table(intermittent=c(rep(0, 300), rep(1, 300)))
#' task = mlr3::TaskClassif$new(id = "in_dt", backend = in_dt, target ='intermittent')
#' get_oversamp_ratio(task)
#' }
#'
#' @export
get_oversamp_ratio <- function(in_task) {
  return(
    in_task$data()[, .N, by=get(in_task$target_names)] %>%
      setorder(N) %>%
      .[, list(minoclass=get[1], ratio=N[2]/N[1])]
  )
}

## create baselearners -----
create_baselearners <- function(in_task, ncores = parallel::detectCores()-3){
  lrns <- list()
  
  if (is.list(in_task)) {
    in_task <- in_task[[1]]
  }
  if (inherits(in_task, "TaskClassif")) {
    
    #Compute ratio of intermittent to perennial observations
    imbalance_ratio <- get_oversamp_ratio(in_task)$ratio
    
    lrns[['lrn_ranger']] <- mlr3::lrn('classif.ranger',
                                      num.trees = 800,
                                      sample.fraction = 0.632,
                                      replace = FALSE,
                                      splitrule = 'gini',
                                      predict_type = 'prob',
                                      importance = 'impurity_corrected',
                                      respect.unordered.factors = 'order')
    
    set_threads(lrns[["lrn_ranger"]], n = ncores)
    
    po_oversampled <- mlr3pipelines::po("classbalancing", id = "oversampled", adjust = "minor",
                                        reference = "minor", shuffle = TRUE,
                                        ratio = imbalance_ratio)
    
    #Create mlr3 pipe operator to put a higher class weight on minority class
    # po_classweights <- mlr3pipelines::po("classweights", minor_weight = imbalance_ratio)
    #Create graph learners so that oversampling happens systematically upstream of all training
    lrns[['lrn_ranger_oversampled']] <- mlr3pipelines::GraphLearner$new(
      po_oversampled %>>% lrns[['lrn_ranger']])
    # lrns[['lrn_ranger_classweights']] <- mlr3pipelines::GraphLearner$new(
    #   po_classweights %>>% lrns[['lrn_ranger']])
  }
  lrns_out <- lrns[['lrn_ranger_oversampled']]
  return(lrns_out)
}

## set tuning parameters ------
set_tuning <- function(in_learner, in_measure, nfeatures,
                       insamp_nfolds, insamp_neval, insamp_nbatch) {
  
  if (is.list(in_learner)) {
    in_learner <- in_learner[[1]]
  }
  
  #Define paramet space to explore
  regex_tuneset <- function(in_learner) {
    prmset <- names(in_learner$param_set$tags)
    tune_rf <- ParamSet$new(list(
      ParamInt$new(grep(".*mtry", prmset, value=T)[1],
                   lower = floor(0.2*nfeatures),
                   upper = floor(0.8*nfeatures)), #80%  of features
      ParamDbl$new(grep(".*fraction", prmset, value=T),
                   lower = 0.2,
                   upper = 0.8)
    ))
    
    in_split =in_learner$param_set$get_values()[
      grep(".*split(rule|stat)", prmset, value=T)]
    
    if (in_split == 'maxstat') {
      tune_rf$add(
        ParamDbl$new(prmset[grep(".*alpha", prmset)],
                     lower = 0.01, upper = 0.1)
      )
      
    } else if (any(grepl(".*min.node.size", prmset))) {
      tune_rf$add(
        ParamInt$new(prmset[grep(".*min.node.size", prmset)],
                     lower = 1, upper = 10)
      )
    } else if (any(grepl(".*splitstat", prmset))) {
      tune_rf$add(
        ParamDbl$new(prmset[grep(".*alpha", prmset)],
                     lower = 0.01, upper = 0.1)
      )
    }
  }
  
  #Define inner resampling strategy
  rcv_rf = rsmp("cv", folds=insamp_nfolds) #aspatial CV repeated 10 times
  
  #Define termination rule
  evalsn = mlr3tuning::trm("evals", n_evals = insamp_neval) #termine tuning after insamp_neval rounds
  
  if (in_learner$task_type == 'classif') {
    if (inherits(in_measure, 'list')) {
      in_measure <- in_measure$classif
    }
    
    if (grepl('classif[.]cforest$', in_learner$id)) {
      learnertune <- in_learner
    } else if (grepl('classif[.]ranger$', in_learner$id)) {
      learnertune <- AutoTuner$new(learner= in_learner,
                                   resampling = rcv_rf,
                                   measure = in_measure,
                                   search_space = regex_tuneset(in_learner),
                                   terminator = evalsn,
                                   tuner =  tnr("random_search",
                                                batch_size = insamp_nbatch)) #batch_size determines level of parallelism
    } else{
      stop('The classification learner provided is not configurable with this workflow yet...')
    }
  } else if (in_learner$task_type == 'regr') {
    if (inherits(in_measure, 'list')) {
      in_measure <- in_measure$regr
    }
    
    learnertune <- AutoTuner$new(learner= in_learner,
                                 resampling = rcv_rf,
                                 measure = in_measure,
                                 search_space = regex_tuneset(in_learner),
                                 terminator = evalsn,
                                 tuner =  tnr("random_search",
                                              batch_size = insamp_nbatch))
  }
  
  #learnertune$store_tuning_instance = FALSE
  learnertune$id <- in_learner$id
  
  return(learnertune)
}

## set_cvresampling -----
set_cvresampling <- function(rsmp_id, in_task, outsamp_nrep, outsamp_nfolds) {
  #repeated_cv or repeated-spcv-coords
  outer_resampling = rsmp(rsmp_id,
                          repeats = outsamp_nrep,
                          folds = outsamp_nfolds)
  outer_resampling$instantiate(in_task)
  
  return(outer_resampling)
}

# dynamic_resample--------
# subsampling the data for dealing with imbalanced data
dynamic_resample <- function(in_task, in_learner, in_resampling, type,
                             store_models = FALSE) {
  if (is.list(in_learner)) {
    in_learner <- in_learner[[1]]
  }
  
  if (is.list(in_task)) {
    in_task <- in_task[[1]]
  }
  
  if (inherits(in_learner, 'BenchmarkResult')) {
    print(('BenchmarkResults was provided, getting the learner...'))
    in_learner <- in_learner$learners$learner[[1]]
  }
  
  #Make sure autotuner matches task (adjust mtry)
  if (inherits(in_learner, 'AutoTuner')) {
    in_learner <- reset_tuning(in_autotuner = in_learner,
                               in_task = in_task)
  }
  
  if ((in_learner$task_type == 'classif' & type=='classif') |
      (in_learner$task_type == 'regr' & type=='regr')) {
    resmp_rs <- mlr3::resample(learner = in_learner, task = in_task,
                               resampling = in_resampling, store_models = store_models
    )
    return(resmp_rs)
  }
}
# weighted_sd -----
weighted_sd <- function(x, w=NULL, na.rm=FALSE) {
  if (na.rm) {
    x <-  na.omit(x)
    if (length(w) > length(x)) {
      w <- w[-which(is.na(x))]
    }
  }
  
  if (length(w)==0) {
    w <- rep(1, length(x))
  }
  
  #Compute weighted standard deviation
  return(sqrt(sum((w) * (x - weighted.mean(x, w)) ^ 2) / (sum(w) - 1)))
}

# combine_bm -----------
combine_bm <- function(in_resampleresults, write_qs = NULL, inp_resdir = NULL) {
  #When tried as_benchmark_result.ResampleResult, got "Error in setcolorder(data, slots) :
  # x has some duplicated column name(s): uhash. Please remove or rename the
  # duplicate(s) and try again.". SO use this instead
  if (!file.exists(inp_resdir)) {
    dir.create(inp_resdir, recursive = TRUE)
  }
  
  print('Converting to benchmark results...')
  if (length(in_resampleresults) > 1) {
    bmres_list <- lapply(
      in_resampleresults[!sapply(in_resampleresults, is.null)],
      function(rsmpres) {
        print(rsmpres)
        if (!is.null(rsmpres)) {
          as_benchmark_result(rsmpres)
        }
      })
    #BenchmarkResult$new(rsmpres$data)})
    
    print('Combining...')
    bmrbase = bmres_list[[1]]
    for (i in 2:length(bmres_list)) {
      if (in_resampleresults[[i]]$task$task_type ==
          in_resampleresults[[1]]$task$task_type) {
        print(i)
        bmrbase$combine(bmres_list[[i]])
      } else {
        warning('ResampleResult #', i,
                'is not of the same task type as the first ResampleResult you provided, skipping...')
      }
    }
  } 
  else {
    warning('You provided only one resample result to combine_bm,
            simply returning output from as_benchmark_result...')
    bmrbase = as_benchmark_result(in_resampleresults[[1]])
  }
  print('Done combining, now writing to qs...')
  if (write_qs) {
    out_filen <- paste0('combine_bm', format(Sys.time(), '%Y%m%d%H%M%s'), '.qs')
    out_qs <- file.path(inp_resdir, out_filen)
  }
  
  qs::qsave(bmrbase, out_qs)
  
  return(out_filen)
}
# select_features ------
select_features <- function(in_bm, in_lrnid, in_task, pcutoff, inp_resdir = NULL) {
  
  #If path, read qs
  if (inherits(in_bm, "character")) {
    in_bm <- qs::qread(file.path(inp_resdir, in_bm))
  }
  
  #get desired resampled_results/learner
  if (inherits(in_bm, "BenchmarkResult")) {
    in_rf <- in_bm$filter(learner_ids = in_lrnid)
  } else {
    in_rf <- as_benchmark_result(in_bm)
  }
  
  #Apply feature/variable selection
  vimp <- weighted_vimportance_nestedrf(
    rfresamp = in_rf$resample_result(uhash=unique(as.data.table(in_rf)$uhash)),
    pvalue = TRUE) %>%
    .[,imp_wmeanper := imp_wmean/sum(imp_wmean)]
  
  task_featsel <- in_task$clone()$select(
    vimp[imp_pvalue <= pcutoff, as.character(varnames)])
  task_featsel$id <- paste0(in_task$id, '_featsel')
  
  return(list(in_task, task_featsel))
}
# weighted_vimportance_nestedrf ----------
weighted_vimportance_nestedrf <- function(rfresamp,
                                          pvalue = TRUE, pvalue_permutn = 50) {
  varnames <- rfresamp$task$feature_names
  rfresampdt <- as.data.table(rfresamp)
  
  vimportance_all <- rfresampdt[, extract_impperf_nestedrf(
    in_rflearner = learner, #Extract vimp and perf for each resampling instance
    in_task = task,
    imp=T, perf=T, pvalue=pvalue, pvalue_permutn), by=iteration] %>%
    cbind(., varnames)
  
  ####!!!!!!!!!!Adapt to allow for other measure than classif.bacc!!!!!!!!######
  out_vimportance <- vimportance_all[
    , list(imp_wmean = weighted.mean(importance, classif.bacc), #Compute weighted mean
           imp_wsd =  weighted_sd(x=importance, w=classif.bacc)), #Compute weighted sd
    by=varnames]
  
  if (pvalue) {
    out_vimportance <- cbind(
      out_vimportance,
      vimportance_all[,
                      list(imp_pvalue = weighted.mean(pvalue, classif.bacc)), #Compute weighted mean of pvalue
                      by=varnames][, !'varnames']
    )
  }
  
  return(out_vimportance)
}

# extract_impperf_nestedrf -----
extract_impperf_nestedrf <- function(in_rflearner, in_task,
                                     imp = TRUE, perf = TRUE,
                                     pvalue = TRUE, pvalue_permutn = 50) {
  
  in_task <- in_task[[1]]
  in_rflearner <- in_rflearner[[1]]
  
  if (inherits(in_rflearner, "AutoTuner")) {
    sublrn <- in_rflearner$model$learner
  } else {
    sublrn <- in_rflearner
  }
  
  print(paste0("Computing variable importance for resampling instance hash #",
               sublrn$hash))
  
  outobj <- cbind(
    if (imp) {
      ####################### IF GraphLearner ####################################
      if (inherits(sublrn, "GraphLearner")) {
        
        if ('classif.ranger' %in% names(sublrn$model)) {
          
          if (pvalue == TRUE) {
            in_formula <- as.formula(paste0(in_task$target_names, '~.'))
            
            importance_pvalues(
              sublrn$model$classif.ranger$model,
              method = "altmann",
              num.permutations = pvalue_permutn,
              data = in_task$data(),
              formula= in_formula
            )
            
          } else {
            data.table(importance=sublrn$model$classif.ranger$model$variable.importance)
          }
        }
        
        else if ('classif.cforest' %in% names(sublrn$model)) {
          if (pvalue == TRUE) {
            warning("p_value calculation is only available for ranger classification rf, ignoring p_value.
                    In addition, default parameters were used in partykit::varimp, adjust as needed.")
          }
          data.table(importance=
                       partykit::varimp(sublrn$model$classif.cforest$model,
                                        nperm = 1,
                                        OOB = TRUE,
                                        risk = "misclassification",
                                        conditional = FALSE,
                                        threshold = .2))
          }
      } else { ####################### IF direct model ####################################
        if (pvalue == TRUE) { #If want pvalue associated with predictor variables
          in_formula <- as.formula(paste0(in_task$target_names, '~.'))
          
          importance_pvalues(
            in_rflearner$model,
            method = "altmann",
            num.permutations = pvalue_permutn,
            data = in_task$data(),
            formula= in_formula
          )
        } else { #If pvalue == FALSE
          data.table(importance= in_rflearner$model$learner$importance())
        }
        
      }
    },
    if (perf) {
      perf_id <- in_rflearner$instance_args$measure$id
      outperf <- in_rflearner$tuning_result[, get(perf_id)]
      data.table(outperf) %>% setnames(perf_id)
    }
  )
  
  return(outobj)
}

# dynamic_resamplebm ----
dynamic_resamplebm <- function(in_task, in_bm, in_lrnid, in_resampling, type,
                               inp_resdir = NULL, store_models = FALSE) {
  #If path, read qs
  if (inherits(in_bm, "character")) {
    in_bm <- qs::qread(file.path(inp_resdir, in_bm))
  }
  
  #get desired resampled_results/learner
  in_rf <- in_bm$filter(learner_ids = in_lrnid)
  
  rsmp_out <- dynamic_resample(in_task = in_task,
                               in_learner = in_rf,
                               in_resampling = in_resampling,
                               type = type,
                               store_models = store_models)
  
  return(rsmp_out)
}
# reset_tuning -----
reset_tuning <- function(in_autotuner, in_task, in_lrnid = NULL) {
  if (inherits(in_autotuner, 'list') & !is.null(in_lrnid)) {
    in_autotuner <- in_autotuner[[
      which(unlist(lapply(in_autotuner, function(lrn) {lrn$id == in_lrnid})))
      ]]
  }
  
  tuneargs_ini <- in_autotuner$instance_args
  
  autotuner_new <- set_tuning(in_learner = tuneargs_ini$learner,
                              in_measure = tuneargs_ini$measure,
                              nfeatures = length(in_task$feature_names),
                              insamp_nfolds= tuneargs_ini$resampling$param_set$values$folds,
                              insamp_neval= tuneargs_ini$terminator$param_set$values$n_evals,
                              insamp_nbatch= in_autotuner$tuner$param_set$values$batch_size
  )
  
  return(autotuner_new)
}

# analyze_benchmark -----
analyze_benchmark <- function(in_bm, in_measure, inp_resdir=NULL) {
  
  #If path, read qs
  if (inherits(in_bm, "character")) {
    in_bm <- qs::qread(file.path(inp_resdir, in_bm))
  }
  
  print(paste('It took',
              in_bm$aggregate(mlr3::msr('time_both'))$time_both,
              'seconds to train and predict with the',
              in_bm$aggregate(msr('time_both'))$learner_id,
              'model...'))
  
  bmdt <- as.data.table(in_bm)
  
  if (in_bm$task_type == 'regr') {
    print(in_bm$aggregate(in_measure$regr))
    boxcomp <- mlr3viz::autoplot(in_bm, measure = in_measure$regr)
    
    preds <- lapply(seq_len(bmdt[,.N]), function(rsmp_i) {
      preds <- bmdt$prediction[[rsmp_i]] %>%
        as.data.table %>%
        .[, `:=`(outf = bmdt$iteration[[rsmp_i]],
                 task = bmdt$task[[rsmp_i]]$id,
                 task_type = in_bm$task_type,
                 learner = bmdt$learner[[rsmp_i]]$id)]
      return(preds)
    }) %>%
      do.call(rbind, .)
    
    if (!('prob.1' %in% names(preds)) & 'response' %in% names(preds)) {
      preds[, prob.1 := response]
    }
  }
  
  
  if (in_bm$task_type == 'classif_st') {
    print(in_bm$aggregate(measures=in_measure$classif))
    boxcomp <- mlr3viz::autoplot(in_bm, measure = in_measure$classif)
    
    preds <- lapply(seq_len(bmdt[,.N]), function(rsmp_i) {
      preds <- data.table(outf = bmdt$iteration[[rsmp_i]],
                          task = bmdt$task[[rsmp_i]]$id,
                          learner = bmdt$learner[[rsmp_i]]$id,
                          task_type = in_bm$task_type,
                          pred = list(bmdt$prediction[[rsmp_i]]))
      return(preds)
    }) %>%
      do.call(rbind, .)
  }
  
  tasklearner_unique <- preds[, expand.grid(unique(task), unique(learner))] %>%
    `colnames<-`(c('task', 'learner')) %>%
    setDT
  
  tasklearner_unique[, learner_format := dplyr::case_when(
    learner == 'classif.ranger'~'default RF',
    learner == 'oversample.classif.ranger'~'default RF - oversampled',
    learner == 'classweights.classif.ranger'~'default RF - weighted classes'
  )]
  
  glist <- lapply(1:nrow(tasklearner_unique), function(tsklrn) {
    print(tasklearner_unique[tsklrn,])
    subpred <- preds[task ==tasklearner_unique$task[tsklrn] &
                       learner == tasklearner_unique$learner[tsklrn],]
    
    ggmisclass_out <- ggmisclass_single(in_predictions = subpred)
    
    gout <- ggmisclass_out$plot +
      ggtitle(paste(tasklearner_unique$task[tsklrn],
                    tasklearner_unique$learner_format[tsklrn])) +
      labs(x='Threshold', y='Value')
    
    if (tsklrn < nrow(tasklearner_unique)) {
      gout <- gout +
        theme(legend.position = 'none')
    }
    
    return(list(plot = ggplotGrob(gout),
                interthres_dt = data.table(
                  learner = as.character(tasklearner_unique$learner[tsklrn]),
                  thresh = ggmisclass_out$interthresh
                )
    )
    )
  }) %>%
    unlist(recursive=F)
  
  
  return(list(
    bm_misclasscomp=do.call("grid.arrange",
                            list(grobs=glist[seq(1, length(glist), 2)])), #Get all plots out of the nested list
    bm_boxcomp = boxcomp,
    interthresh_dt = rbindlist(glist[seq(2, length(glist), 2)]) #Get all threshold data.table rows out of nested list
  ))
}

# get_outerrsmp ---------
get_outerrsmp <- function(in_rftuned, spatial_rsp=FALSE) {
  #Adapt whether return resample result or output from selecttrain_rf
  if (inherits(in_rftuned, 'list')) {
    #If there is more than one type of outer resampling
    if (length(in_rftuned$rf_outer$uhash) > 1) {
      #Check which resampling is spatial — only works if one is spatial
      sp_i <- which(unlist(lapply(in_rftuned$rf_outer, function(x) {
        grepl('.*Resampling.*Sp.*', x$resampling$format())
      })))
      
      #If user request that spatial resampling be used
      if (spatial_rsp==TRUE) {
        
        if (length(sp_i)>0) {
          rsmp_res <- in_rftuned$rf_outer[[min(sp_i)]]
        } else { #But if there is no spatial resampling provided
          stop("spatial_rsp==TRUE but the in_rftuned does not include
               any Spatial Resampling")
        }
        #If user didn't request spatial resampling to be used, grab the first
        #resampling that is not spatial
        } else {
          rsmp_res <- in_rftuned$rf_outer[[
            min((1:length(in_rftuned$rf_outer))[-sp_i])]]
        }
      
      #If there is only one type of outer resampling
    } else {
      print("Only one resampling result, ignoring spatial_rsp argument...")
      rsmp_res <- in_rftuned$rf_outer
    }
    #If in_rftuned is already a ResampleResult, simply return it
  } else if (inherits(in_rftuned, "ResampleResult")) {
    rsmp_res <- in_rftuned
  }
  return(rsmp_res)
}

# ggmisclass_single ---------
ggmisclass_single <- function(in_predictions=NULL, in_rftuned=NULL, spatial_rsp=FALSE) {
  #Get predicted probabilities of intermittency for each gauge
  # in_gaugestats[!is.na(cly_pc_cav), intermittent_predprob :=
  #                 as.data.table(in_predictions)[order(row_id), mean(prob.1), by=row_id]$V1]
  #Get misclassification error, sensitivity, and specificity for different classification thresholds
  #i.e. binary predictive assignment of gauges to either perennial or intermittent class
  
  #If provided resampling results rather than prediction table, extract 
  if (!is.null(in_rftuned)) {
    rsmp_res <- get_outerrsmp(in_rftuned, spatial_rsp=spatial_rsp)
    in_predictions <- rsmp_res$prediction()
  }
  
  #Get confusion matrices for range of thresholds (i.e., probability of flow intermittence
  #above which a watercourse is classified as non-perennial)
  threshold_confu_dt <- ldply(seq(0,1,0.01), threshold_misclass, in_predictions) %>%
    setDT
  
  #Get classification threshold at which sensitivity and specificity are the most similar
  balanced_thresh <- threshold_confu_dt[which.min(abs(spec-sens)),]
  print(paste('Sensitivity =', round(balanced_thresh$sens,2),
              'and Specificity =', round(balanced_thresh$spec,2),
              'at a classification threshold of', balanced_thresh$i))
  
  #Plot trends in confusion matrix metrics with increasing threshold
  gout <- ggplot(melt(threshold_confu_dt, id.vars='i'),
                 aes(x=i, y=value, color=variable, linetype=variable)) +
    geom_line(size=1.2) +
    geom_vline(xintercept=balanced_thresh$i, alpha=1/2) +
    geom_hline(yintercept=balanced_thresh$spec, alpha=1/2) +
    annotate('text', x=(balanced_thresh$i), y=0.4,
             label=balanced_thresh$i, angle=-90) +
    annotate('text', x=0.9, y=(balanced_thresh$spec),
             label=round(balanced_thresh$sens,2)) +
    scale_x_continuous(expand=c(0,0), name='Threshold') +
    scale_y_continuous(expand=c(0,0), name='Value') +
    scale_color_brewer(palette='Dark2',  #colorblind friendly
                       labels=c('Misclassification rate',
                                'Sensitivity (true positives)',
                                'Specificity (true negatives)')) +
    theme_bw()
  
  #Plot it
  return(list(plot = gout,
              interthresh = balanced_thresh$i))
}

# threshold_misclass ------
threshold_misclass <- function(i=0.5, in_preds) {
  #---- Get confusion matrix ----
  if (inherits(in_preds, 'PredictionClassif')) {
    confu <- as.data.table(in_preds$set_threshold(1-i)$confusion) #Get confusion matrix directly
  }
  
  if (is.data.table(in_preds)) {
    #If task associated with predictions is a classification and has records
    if (in_preds[task_type == 'classif_st',.N] > 0) {
      #For each CV repetition:
      #   1. set the probability threshold to compute a confusion matrix to 1-i
      #     (i being the threshold to classify something as 1, set_threshold
      #     being based on prob.0, not prob.1)
      #   2. Compute confusion matrix
      confu <- in_preds[, as.data.table(pred[[1]]$set_threshold(1-i)$confusion),
                        by=outf] %>%
        #Aggregate confusion matrices across repetitions
        .[, .(N=sum(N)), by=.(response, truth)]
    }
    
    #If task associated with predictions is a regression and has records
    if (in_preds[task_type == 'regr', .N] > 0) {
      #Reclassify continuous predictions into binary response across all records
      confu <- in_preds[, response := fifelse(prob.1>=i, '1', '0')] %>%
        .[, truth := as.character(truth)] %>%
        #Create aggregate confusion matrix
        .[, .N, by=.(response, truth)]
    }
  }
  
  #---- Compute statistics based on confusion matrix and format into data.table----
  outvec <- data.table(
    i,
    misclas = confu[truth != response, sum(N)] / confu[, sum(N)],
    sens = confu[truth == '1' & response == '1', N] / confu[truth=='1', sum(N)],
    spec  = confu[truth=='0' & response==0, N]/confu[truth=='0', sum(N)]
  )
  return(outvec)
}

# selecttrain_rf --------
selecttrain_rf <- function(in_rf, in_learnerid=NULL, in_task = NULL,
                           insamp_nfolds =  NULL, insamp_nevals = NULL) {
  
  outlist <- list()
  
  ######### Prepare autotuner for full training ####################
  # If a ResampleResult was provided
  if (inherits(in_rf, 'ResampleResult')) {
    in_bmsel <- in_rf$clone()
    iter_selected <- in_bmsel$score() %>% 
      as.data.table() %>%
      .[order(classif.ce)] %>%
      .[1,iteration]
    lrn_autotuner <- in_bmsel$learners[[iter_selected]]
    in_task <- in_bmsel$task
    outlist[['rf_outer']] <- in_rf
    
    # If a BenchmarkResult was provided
  } else if (inherits(in_rf, 'BenchmarkResult')) {
    in_bmsel <- in_rf$clone()$filter(learner_ids = in_learnerid,
                                     task_id = in_task)
    
    lrn_autotuner <- in_bmsel$clone()$learners$learner[[1]]
    in_task <-in_bmsel$tasks$task[[1]]
    
    #Return outer sampling object for selected model (or list of outer sampling objects)
    uhashes <- unique(as.data.table(in_bmsel)$uhash)
    if (length(uhashes) == 1) {
      outlist[['rf_outer']] <- in_bmsel$resample_result(uhash=uhashes)
    } else {
      outlist[['rf_outer']] <- lapply(uhashes, function(x) {
        in_bmsel$resample_result(uhash=x)
      })
    }
  } else if (inherits(in_rf, 'AutoTuner')) {
    lrn_autotuner <- in_rf
  }
  
  if (!is.null(insamp_nfolds)) {
    lrn_autotuner$instance_args$resampling$param_set$values$folds <- insamp_nfolds
  }
  
  if (!is.null(insamp_nevals)) {
    lrn_autotuner$instance_args$terminator$param_set$values$n_evals <- insamp_nevals
  }
  
  ######### Train it ####################
  # lrn_autotuner$model$learner$param_set$values = mlr3misc::insert_named(
  #   lrn_autotuner$model$learner$param_set$values,
  #   list(classif.ranger.importance = 'permutation')
  # )
  lrn_autotuner$train(in_task)
  
  outlist[['task']] <- in_task
  outlist[['rf_inner']] <- lrn_autotuner
  
  return(outlist)
}

# ----------------Functions for step 2 RF ------------

create_task_steptwo <- function(in_tbl, breaks, labels){
  colnames(in_tbl) <- make.names(names(in_tbl),
                                 unique = TRUE)
  
  in_tbl = in_tbl[target>0,][,c('gaugeid', 'dates') := NULL]
  
  target_class <- classify_target(in_tbl,
                                  target_col = 'target',
                                  breaks = breaks, 
                                  labels = labels)
  in_tbl_editted <- in_tbl[, target_class := target_class][,target :=NULL]
  setnames(in_tbl_editted, 'target_class', 'target')
  
  if (!inherits(in_tbl[,target], "factor")) {
    in_tbl[, target := as.factor(target)]
  }
  
  task <- mlr3spatiotempcv::TaskClassifST$new(
    id = "multi_classes",
    backend = in_tbl_editted,
    target = "target",
    coordinate_names = c("X", "Y"))
  
  
  return(task)
}

create_baselearners_steptwo <- function(in_task, ncores = parallel::detectCores()-3){
  lrns <- list()
  
  if (is.list(in_task)) {
    in_task <- in_task[[1]]
  }
  if (inherits(in_task, "TaskClassif")) {
    
    #Compute ratio of intermittent to perennial observations
    # imbalance_ratio <- get_oversamp_ratio(in_task)$ratio
    
    lrns[['lrn_ranger']] <- mlr3::lrn('classif.ranger',
                                      num.trees = 800,
                                      sample.fraction = 0.632,
                                      replace = FALSE,
                                      splitrule = 'gini',
                                      predict_type = 'prob',
                                      importance = 'impurity_corrected',
                                      respect.unordered.factors = 'order')
    
    set_threads(lrns[["lrn_ranger"]], n = ncores)
    
    po_oversample_major <- po("classbalancing", id='oversample_major', reference = 'major',
                              shuffle = FALSE, ratio = 1)
    
    # po_smote <-
    #   po("colapply", id = "int_to_num",
    #      applicator = as.numeric, affect_columns = selector_type("integer")) %>>%
    #   po("smote", dup_size = 2) %>>%
    #   po("colapply", id = "num_to_int",
    #      applicator = function(x) as.integer(round(x, 0L)), affect_columns = selector_type("numeric"))
    
    # po_smote <- po("smote", id='over_smote')
    
    lrns[['lrn_ranger_oversampled']] <- mlr3pipelines::GraphLearner$new(
      po_oversample_major %>>% lrns[['lrn_ranger']])
    
    # lrns[['lrn_ranger_smote']] <- mlr3pipelines::GraphLearner$new(
    #   po_smote %>>% lrns[['lrn_ranger']])
    
  }
  lrns_out <- lrns[['lrn_ranger_oversampled']]
  
  return(lrns)
}

select_features_step2 <- function(in_bm, in_taskid, in_task, pcutoff, inp_resdir = NULL) {
  
  #If path, read qs
  if (inherits(in_bm, "character")) {
    in_bm <- qs::qread(file.path(inp_resdir, in_bm))
  }
  
  #get desired resampled_results/learner
  if (inherits(in_bm, "BenchmarkResult")) {
    in_rf <- in_bm$filter(learner_ids = 'classif.ranger', task_ids = in_taskid)
  } else {
    in_rf <- as_benchmark_result(in_bm)
  }
  
  #Apply feature/variable selection
  vimp <- weighted_vimportance_nestedrf(
    rfresamp = in_rf$resample_result(uhash=unique(as.data.table(in_rf)$uhash)),
    pvalue = TRUE) %>%
    .[,imp_wmeanper := imp_wmean/sum(imp_wmean)]
  
  task_featsel <- in_task$clone()$select(
    vimp[imp_pvalue <= pcutoff, as.character(varnames)])
  task_featsel$id <- paste0(in_task$id, '_featsel')
  
  return(list(in_task, task_featsel))
}

set_cvresampling_step2 <- function(rsmp_id, in_task, outsamp_nrep, outsamp_nfolds) {
  
  if (is.list(in_task)) {
    in_task <- in_task[[1]]
  }
  #repeated_cv or repeated-spcv-coords
  outer_resampling = rsmp(rsmp_id,
                          repeats = outsamp_nrep,
                          folds = outsamp_nfolds)
  # outer_resampling = Map(function(task) outer_resampling$instantiate(task),
  #                        task = in_task)
  outer_resampling$instantiate(in_task)
  
  return(outer_resampling)
}
