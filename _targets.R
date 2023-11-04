library(targets)
library(tarchetypes)
tar_option_set(packages = c("sf", "data.table", "tidyverse", "here"),
               memory = "transient", garbage_collection = TRUE, format = "qs")

tar_source("src/used_libraries.R")
tar_source("src/pre_processing_functions.R")
tar_source("src/modeling_functions.R")
tar_source("src/producing_figures.R")


# define the data path 
rootpath = here::here()
datapath = file.path(rootpath, "Data")
shp_path = file.path(rootpath, "Data/shp/eu_stations/european_gaugingstations.shp")
resultspath = file.path(rootpath, "results")
check_plot_path = file.path(resultspath, 'ts_plots')
preprocessed_path = file.path(resultspath, "preprocessed")
check_plot_set1 = file.path(resultspath, 'ts_plots/set1')
check_plot_set2 = file.path(resultspath, 'ts_plots/set2')
shapefiles_path = file.path(resultspath, 'shapefiles')
model_path = file.path(resultspath, 'modeldir')
reachoutput_path = file.path(resultspath, 'reach_output')

lapply(list(resultspath, check_plot_path, preprocessed_path,
            check_plot_set1, check_plot_set2, shapefiles_path, model_path,
            reachoutput_path), function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
})

# Pathes used when the large rasters are available in the workflow 
watergap_raster_HR_path = file.path("/home/home1/gm/projects/DRYvER/03_data/13_predictors", 
                                    "Downscaled_WaterGAP_raster")

path_LR = "/home/home1/gm/projects/DRYvER/03_data/14_Mathis_data/LRpredictors"
reachesdatpath = '/home/home1/gm/projects/DRYvER/04_code/04_statistical_modeling/Apply_nets_eu/data'
# define start and end dates
start_date = as.Date("1981-01-01")
end_date = as.Date("2019-12-31")
# ----------------------- pre-processing: select interested European stations -----------------
plan_preprocess = tar_plan(
  tar_target(
    name = 'spanish_rawfiles',
    download_spanish_stations(file.path(rootpath, "Data/targets/raw", "Spain_SAIH"))
  ),
  
  tar_target(
    name = "input_files",
    command = list(
      list.files(
        path = file.path(datapath, "targets/raw/2023-04-20_08-00_Europe1"),
        pattern = "txt",
        full.names = TRUE
      ),
      list.files(
        path = file.path(datapath, "targets/raw/corsica"),
        pattern = "csv",
        full.names = TRUE
      ),
      file.path(datapath,
                "targets/raw/harmonized_dataset_cleaned_new.csv"
                ),
      list.files(
        path = file.path(datapath, "targets/raw/Smires/uuid"),
        pattern = "csv",
        full.names = TRUE
        ),
      list.files(
        path = file.path(datapath, "targets/raw/Data_IT/emr"),
        pattern = "csv",
        full.names = TRUE
        ),
      list.files(
        path = file.path(datapath, "targets/raw/Data_IT/ispra"),
        pattern = "csv",
        full.names = TRUE
        ),
      list.files(
        path = file.path(datapath, "targets/raw/Data_IT/lom"),
        pattern = "csv",
        full.names = TRUE
        ),
      list.files(
        path = file.path(datapath, "targets/raw/Data_IT/sar"),
        pattern = "csv",
        full.names = TRUE
        )
    )
  ),
  tar_target(
    name = "grdc_stations",
    command = select_new_grdc_stations(path = input_files[[1]],
                                       shp_path = shp_path)
  ),
  tar_target(
    name = "corsica_stations",
    command = select_corsica_stations(path = input_files[[2]],
                                      shp_path = shp_path)
  ),
  tar_target(
    name = "grdc_old_stations",
    command = select_old_grdc_stations(path = input_files[[3]],
                                       shp_path = shp_path,
                                       old_grdc_stations_id = grdc_stations$old_grdc_stations_id)
  ),
  tar_target(
    name = "gsim_stations",
    command = select_gsim_stations(path = input_files[[3]],
                                   shp_path = shp_path)
  ),
  tar_target(
    name = "rbis_stations",
    command = select_gsim_stations(path = input_files[[3]],
                                   shp_path = shp_path,
                                   dataset_name = "RBIS")
  ),

  tar_target(
    name = "smires_stations",
    command = select_smires_stations(path = input_files[[4]],
                                     shp_path = shp_path)
  ),
  tar_target(
    name = "emr_stations",
    command = select_italian_emr_stations(path = input_files[[5]],
                                          shp_path = shp_path)
  ),
  tar_target(
    name = "ispra_stations",
    command =  select_italian_ispra_stations(path = input_files[[6]],
                                             shp_path = shp_path)
    
  ),
  tar_target(
    name = "arpal_stations",
    command = select_italian_arpal_stations(path = input_files[[7]],
                                            shp_path = shp_path)
  ),
  tar_target(
    name = "arpas_stations",
    command =  select_italian_arpas_stations(path = input_files[[8]],
                                             shp_path = shp_path)
  ),
  
  tar_target(
    name = 'spanish_stations',
    command = select_spanish_stations(in_paths = spanish_rawfiles$paths,
                                      in_metadata =  spanish_rawfiles$metadata
    )
  ),
  # tar_target(
  #       name = "write_files",
  #       command = {
  #         data.table::fwrite(
  #           corsica_stations$output_ts, file = file.path(datapath,"preprocessed/corsica_stations.csv"))
  #         data.table::fwrite(
  #           grdc_stations$output_dt, file = file.path(datapath,"preprocessed/grdc_updated_stations.csv"))
  #         data.table::fwrite(
  #           grdc_old_stations, file = file.path(datapath, "preprocessed/grdc_old_139stations.csv"))
  #         data.table::fwrite(
  #           gsim_stations, file = file.path(datapath, "preprocessed/gsim_stations.csv"))
  #         data.table::fwrite(
  #           rbis_stations, file = file.path(datapath, "preprocessed/rbis_stations.csv"))
  #         data.table::fwrite(
  #           smires_stations$output_ts, file = file.path(datapath, "preprocessed/smires_stations.csv"))
  #         data.table::fwrite(
  #           emr_stations$output_ts, file = file.path(datapath, "preprocessed/emr_stations.csv"))
  #         data.table::fwrite(
  #           ispra_stations$output_ts, file = file.path(datapath,"preprocessed/ispra_stations.csv"))
  #         data.table::fwrite(
  #           arpal_stations$output_ts, file = file.path(datapath, "preprocessed/arpal_stations.csv"))
  #         data.table::fwrite(
  #           arpas_stations$output_ts, file = file.path(datapath, "preprocessed/arpas_stations.csv"))
  #       }, format = "file"
  # ),
  # 
  # tar_target(
  #   ts_plots,
  #   lapply(list(corsica_stations$output_ts,
  #               grdc_stations$output_dt,
  #               grdc_old_stations$output_dt,
  #               smires_stations$output_ts,
  #               emr_stations$output_ts,
  #               ispra_stations$output_ts,
  #               arpal_stations$output_ts,
  #               arpas_stations$output_ts,
  #               spanish_stations$output_ts,
  #               gsim_stations$output_ts,
  #               rbis_stations$output_ts),
  #   function(in_dt) {
  #     lapply(unique(in_dt$gaugeid),
  #            function(id) plot_daily_q_timeseries(in_dt[gaugeid == id,],
  #                                                 outdir=check_plot_path))
  #   }
  #   )
  # ),
  # tar_target(
  #   name = 'split_plots',
  #   command = split_file_random(source_path = check_plot_path,
  #                               set1_path = check_plot_set1,
  #                               set2_path = check_plot_set2)
  # )
  tar_target(
    name = final_targets,
    command = {
      list(corsica_stations$output_ts,
           grdc_stations$output_dt,
           grdc_old_stations$output_ts,
           smires_stations$output_ts,
           emr_stations$output_ts,
           ispra_stations$output_ts,
           arpal_stations$output_ts,
           arpas_stations$output_ts,
           spanish_stations$output_ts,
           gsim_stations$output_ts,
           rbis_stations$output_ts) %>% 
        remove_records(.)
    }
  ),
  tar_target(
    name = "final_gauges",
    command = final_gauge_shp(shp_path,
                              gaugeids_removed = final_targets$gaugeids_removed)
  )
)

# ----------------------- pre-processing: computing the predictors -----------------
plan_compute_predictors = tar_plan(
  
  # extract the HR, LR, and static predictors from large files that cannot be stores
  # in clouds, so we decided to comment this section in the workflow. on the other hand,
  # the predictors for the stations are stored in csv format.
  
  # tar_target(
  #   name = "watergap_hr_ts",
  #   command = extract_HR_streamflow(path = watergap_raster_HR_path,
  #                                   gauges_shp = final_gauges,
  #                                   out_dir = file.path(datapath, "predictors/HR"))
  # ),
  # tar_target(
  #   name = "extract_qrdif_ql_ratio",
  #   command = extract_LR_predictors(path = file.path(path_LR, "qrdif_ql_ratio_mon_raster"),
  #                                   gauges_shp = final_gauges,
  #                                   variable_name = "qrdif_ql_ratio",
  #                                   out_dir = file.path(datapath, "predictors/LR"))
  # ),
  # tar_target(
  #   name = "extract_wetdays",
  #   command = extract_LR_predictors(path = file.path(path_LR, "Prec_wetdays_raster"),
  #                                   gauges_shp = final_gauges,
  #                                   variable_name = "wetdays",
  #                                   out_dir = file.path(datapath, "predictors/LR"))
  # ),
  tar_target(
    name = "hr_predictors",
    command = compute_highres_predictors(
      watergap_raw_path = file.path(datapath, "predictors/HR",
                                    "waterGAP_streamflow_stations.csv"))
  ),
  tar_target(
    name = 'static_predictors',
    command = compute_static_predictors(
      in_dt_path = file.path(datapath, "eu_nets", "attri_table_eu_nets.csv"),
      gauge_shp = final_gauges
    )
  ),
  tar_target(
    name = "model_data",
    command = combined_predictors(in_target = final_targets$output_ts,
                                  in_hr_pred = hr_predictors,
                                  path_lr_pred = file.path(datapath, "predictors/LR"),
                                  in_static_pred = static_predictors)
  )
)

# ------------------------- Modeling development- step ONE --------------------------


plan_model_stepone = tar_plan(
  # tar_target(
  #   name = "n",
  #   command = {
  #     n = sample(1:nrow(model_data), 40000)
  #   }
  # ),
  tar_target(
    name = "in_task",
    command = create_task_stepone(in_tbl = model_data[,-c(1, 2)])
  ),
  tar_target(
    name = "measures",
    command = list(classif = msr("classif.bacc"),
                   regr = msr("regr.mae"))
  ),
  tar_target(
    name = "baselearners",
    command = create_baselearners(in_task = in_task, ncores = 18)
  ),
  tar_target(
    name = "seplearners",
    command = tar_read(baselearners),
    iteration = "list"
  ),
  tar_target(
    name = "autotuning",
    command = set_tuning(in_learner = seplearners,
               in_measure = measures,
               nfeatures = length(in_task$feature_names),
               insamp_nfolds = 3, insamp_neval = 15,
               insamp_nbatch = parallel::detectCores(logical = FALSE) - 4),
    pattern = map(seplearners)
  ),
  tar_target(
    name = "resamplingset",
    command = set_cvresampling(rsmp_id = "repeated_cv",
                               in_task = in_task,
                               outsamp_nrep = 2,
                               outsamp_nfolds = 3)
  ),
  tar_target(
    name = "rfresampled_classif",
    command = dynamic_resample(in_task = in_task,
                               in_learner = autotuning,
                               in_resampling = resamplingset,
                               store_models = TRUE,
                               type = "classif"),
    pattern = map(autotuning),
    iteration = "list"

  ),
  tar_target(
    name = "rfbm_classif",
    command = combine_bm(in_resampleresults = rfresampled_classif,
                         write_qs = T,
                         inp_resdir = file.path(resultspath, "store_premodels_stepone"))
  ),
  tar_target(
    name = "selected_learner",
    command = "oversampled.classif.ranger"
  ),
  tar_target(
    name = "tasks_featsel",
    command = select_features(
      in_bm = rfbm_classif,
      in_lrnid =  selected_learner,
      in_task = in_task,
      pcutoff = 0.05,
      inp_resdir = file.path(resultspath, "store_premodels_stepone")
    )
  ),
  tar_map(
    values = tibble(in_strategy = c('repeated_cv', "repeated_spcv_coords"),
                    in_outrep = c(2, 1),
                    in_outfolds = c(3, 10),
                    names = c('cv', 'spcv')),
    names =  names,
    tar_target(
      name = "featsel",
      command = set_cvresampling(rsmp_id = in_strategy,
                                 in_task = in_task,
                                 outsamp_nrep = in_outrep,
                                 outsamp_nfolds = in_outfolds)
  ),
  unlist = FALSE
  ),
  tar_target(
    name = "res_featsel_cv",
    command = dynamic_resamplebm(in_task = tasks_featsel[[2]],
                                 in_bm = rfbm_classif,
                                 in_lrnid =  selected_learner,
                                 in_resampling = featsel_cv,
                                 store_models = TRUE,
                                 inp_resdir = file.path(resultspath, "store_premodels_stepone"),
                                 type = 'classif')
  ),
  tar_target(
    name = "res_featsel_spcv",
    command = dynamic_resamplebm(in_task = tasks_featsel[[2]],
                                 in_bm = rfbm_classif,
                                 in_lrnid =  selected_learner,
                                 in_resampling = featsel_spcv,
                                 store_models =  TRUE,
                                 inp_resdir = file.path(resultspath, "store_premodels_stepone"),
                                 type = 'classif')
  ),
  tar_target(
    name = "rfeval_featsel",
    c(res_featsel_cv, res_featsel_spcv)
  ),
  tar_target(
    name = "rfbm_featsel",
    command = analyze_benchmark(in_bm = rfeval_featsel,
                                in_measure = measures)
  ),
  tar_target(
    name = "rftuned",
    command = selecttrain_rf(in_rf = res_featsel_cv,
                             in_learnerid = selected_learner,
                             in_task = "binary_class")
  ),
  tar_target(
    name = predvars,
    command = predname_df(in_task = rftuned$task)
  ),
  tar_target(
    name = 'vimp_plot',
    command = ggvimp(in_rftuned = rftuned, in_predvars = predvars,
                     varnum=22, spatial_rsp = FALSE)
  ),
   tar_target(
     name = "pd_plot",
     command = ggpartialdep(in_rftuned=rftuned,
                            in_predvars=predvars,
                            colnums=1:27, nvariate=1, nodupli = FALSE,
                            grid_resolution = 20, parallel = FALSE, spatial_rsp = FALSE)
   )
)

# ------------------ Modeling step2 with 4 categories test --------------
plan_model_step2_fourclasses = tar_plan(
  tar_target(
    name = "in_task_step2_fourclass",
    command = create_task_steptwo_test(in_tbl = model_data,
                                       breaks = c(-Inf, 5,15, 29, 32),
                                       labels = c("1", "2", "3", "4"))
  ),
  tar_target(
    name = "measures",
    command = list(classif = msr("classif.bacc"),
                   regr = msr("regr.mae"))
  ),
  tar_target(
    name = "baselearners_step2_fourclass",
    command = create_baselearners_steptwo(in_task = in_task_step2_fourclass,
                                                ncores = 18)
  ),
  tar_target(
    name = "seplearners_step2_fourclass",
    command = tar_read(baselearners_step2_fourclass),
    iteration = "list"
  ),
  tar_target(
    name = "autotuning_step2_fourclass",
    command = set_tuning(in_learner = seplearners_step2_fourclass,
                         in_measure = measures,
                         nfeatures = length(in_task_step2_fourclass$feature_names),
                         insamp_nfolds = 3, insamp_neval = 55,
                         insamp_nbatch = parallel::detectCores(logical = FALSE) - 1),
    pattern = map(seplearners_step2_fourclass)
  ),
  tar_target(
    name = "resamplingset_step2_fourclass",
    command = set_cvresampling_step2(rsmp_id = "repeated_cv",
                                     in_task = in_task_step2_fourclass,
                                     outsamp_nrep = 2,
                                     outsamp_nfolds = 3)
  ),
  tar_target(
    name = "rfresampled_classif_step2_fourclass",
    command = dynamic_resample(in_task = in_task_step2_fourclass,
                               in_learner = autotuning_step2_fourclass,
                               in_resampling = resamplingset_step2_fourclass,
                               store_models = TRUE,
                               type = "classif"),
    pattern = map(autotuning_step2_fourclass),
    iteration = "list"
    
  ),
  tar_target(
    name = "rfbm_classif_step2_fourclass",
    command = combine_bm(in_resampleresults = rfresampled_classif_step2_fourclass,
                         write_qs = T,
                         inp_resdir = file.path(resultspath, "store_premodels_steptwo"))
  ),
  tar_target(
    name = "selected_learner_step2_fourclass",
    command = "oversample_major.classif.ranger"
  ),
  tar_target(
    name = "tasks_featsel_step2_fourclass",
    command = select_features(
      in_bm = rfbm_classif_step2_fourclass,
      in_lrnid =  selected_learner_step2_fourclass,
      in_task = in_task_step2_fourclass,
      pcutoff = 0.05,
      inp_resdir = file.path(resultspath, "store_premodels_steptwo")
    )
  ),
  tar_map(
    values = tibble(in_strategy_step2 = c('repeated_cv', "repeated_spcv_coords"),
                    in_outrep_step2 = c(2, 1),
                    in_outfolds_step2 = c(3, 10),
                    names = c('cv', 'spcv')),
    names =  names,
    tar_target(
      name = "featsel_step2_fourclass",
      command = set_cvresampling(rsmp_id = in_strategy_step2,
                                 in_task = in_task_step2_fourclass,
                                 outsamp_nrep = in_outrep_step2,
                                 outsamp_nfolds = in_outfolds_step2)
    ),
    unlist = FALSE
  ),
  # tar_target(
  #   name = "res_featsel_cv_step2_fourclass",
  #   command = dynamic_resamplebm(in_task = tasks_featsel_step2_fourclass,
  #                                in_bm = rfbm_classif_step2_fourclass,
  #                                in_lrnid =  selected_learner_step2_fourclass,
  #                                in_resampling = featsel_step2_fourclass_cv,
  #                                store_models = TRUE,
  #                                inp_resdir = file.path(resultspath, "store_premodels_stepone"),
  #                                type = 'classif')
  # 
  # )
  # tar_target(
  #   name = "res_featsel_spcv_step2_fourclass",
  #   command = dynamic_resamplebm_step2(in_task = tasks_featsel_step2_fourclass[[2]],
  #                                      in_bm = rfbm_classif_step2_fourclass,
  #                                      in_lrnid =  'classif.ranger',
  #                                      in_taskid = selected_method_step2_fourclass,
  #                                      in_resampling = featsel_step2_fourclass_spcv,
  #                                      store_models =  TRUE,
  #                                      inp_resdir = file.path(resultspath, "store_premodels_steptwo"),
  #                                      type = 'classif')
  # ),
  # tar_target(
  #   name = "rfeval_featsel_step2_fourclass",
  #   c(res_featsel_cv_step2_fourclass, res_featsel_spcv_step2_fourclass)
  # ),
  # # tar_target(
  # #   name = "rfbm_featsel_step2_fourclass",
  # #   command = analyze_benchmark(in_bm = rfeval_featsel_step2_fourclass,
  # #                               in_measure = measures)
  # # ),
  tar_target(
    name = "rftuned_step2_fourclass",
    command = selecttrain_rf(in_rf = rfresampled_classif_step2_fourclass[[2]],
                             in_learnerid = 'classif.ranger',
                             in_task = "multi_class")
  ),
  tar_target(
    name = predvars_step2_fourclass,
    command = predname_df(in_task = rftuned_step2_fourclass$task,
                          in_task_all = in_task_step2_fourclass)
  ),
  tar_target(
    name = 'vimp_plot_step2_fourclass',
    command = ggvimp(in_rftuned = rftuned_step2_fourclass,
                     in_predvars = predvars_step2_fourclass,
                     varnum=23, spatial_rsp = FALSE)
  )
  # tar_target(
  #   name = "pd_plot_Step2_fourclass",
  #   command = ggpartialdep(in_rftuned=rftuned_step2_fourclass,
  #                          in_predvars=predvars_step2_fourclass,
  #                          colnums=1:27, nvariate=1, nodupli = FALSE,
  #                          ngrid = 20, parallel = FALSE, spatial_rsp = FALSE)
  # )
)

# ------------------ Plan for applying models to European reaches ----------
plan_applyingmodels = tar_plan(
  tar_target(
    name = 'eu_reaches_status_dt',
    command = runmodels_over_period(path_model1 = file.path(model_path, "rftuned_step1.qs"),
                                    path_model2 = file.path(model_path, "rftuned_step2.qs"),
                                    path_static = file.path(reachesdatpath, "Statics", "static_preds_net_eu.fst"),
                                    path_LR = file.path(reachesdatpath, "LR"),
                                    path_HR = file.path(reachesdatpath, "HR"),
                                    outdir = reachoutput_path)
  )
)
# ------------------ Plan for producing Figures -------------------
plan_figures = tar_plan(
  tar_target(
    name = 'validate_DS_nselognse', # the format is geospatial and used to create Figure 2 as well
    command = compute_logNSE_NSE_dswatergap(observed_streamflow = final_targets$output_ts,
                                            ds_watergap = hr_predictors,
                                            shp = final_gauges)
  ),
  tar_target(
    name = 'save_validation_figure2',
    command = {
      tar_read(validate_DS_nselognse) %>%
        sf::st_write(., dsn = file.path('results/shapefiles', 'logNSE_NSE_obs_vs_DSwatergap.shp'))
    }
  ),
  tar_target(
    name = 'boxplot_figure3',
    command = ggboxplot_lognse_fig3(in_data = validate_DS_nselognse)
  ),
  tar_target(
    name = 'boxplot_figure5S',
    command = ggboxplot_intermittent_figS5(in_nse_data=validate_DS_nselognse,
                                           in_model_data=model_data)
  ),
  tar_target(
    name = 'boxplot_figure6S',
    command = ggboxplot_perennial_figS6(in_nse_data=validate_DS_nselognse,
                                           in_model_data=model_data)
  ),
  tar_target(
    name = 'ratio_intermittency_step1_figure4ab',
    command = compute_noflow_ratio(in_model=rftuned,
                                   gauges_shp=final_gauges,
                                   outdir='results/shapefiles')
  ),
  tar_target(
    name = 'corr_nse_step1_figure4c',
    command = compute_corr_nse(in_model=rftuned,
                               gauges_shp=final_gauges,
                               outdir='results/shapefiles')
  ),
  tar_target(
    name = 'confmat_step2_figure5',
    command = ggplotConfusionMatrix_fig5(rftuned_step2=rftuned_step2_fourclass,
                                         model_data=model_data)
  ),
  tar_target(
    name = 'varimp_bothstpes_figure7',
    command = ggvimp(in_rftuned_step1=rftuned,
                     in_rftuned_step2=rftuned_step2_fourclass,
                     in_predvars=predvars)
  ),
  tar_target(
    name = 'partialdepplot_figure6S',
    command = ggpartialdep(in_rftuned=rftuned,
                           in_predvars=predvars,
                           model_data = model_data,
                           colnums=1:23, nvariate=1, nodupli = FALSE,
                           ngrid = 40, parallel = FALSE, spatial_rsp = FALSE)
  )
) 

# ------------------ Pipeline of the workflow's plans -------------
list(
  plan_preprocess,
  plan_compute_predictors,
  plan_model_stepone,
  plan_model_step2_fourclasses,
  plan_applyingmodels,
  plan_figures
)
