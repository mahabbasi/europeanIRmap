library(targets)
library(tarchetypes)
tar_option_set(packages = c("sf", "data.table", "tidyverse", "here"),
               memory = "transient", garbage_collection = TRUE, format = "qs",
               future::plan(future::multisession, workers = 10))
tar_source("src/pre_processing_functions.R")
tar_source("src/used_libraries.R")

# define the data path 
rootpath = here::here()
datapath = file.path(rootpath, "data")
resultpath =  file.path(rootpath, "result")
shp_path = file.path(rootpath, "data/shp/eu_stations/spatial_joint_eu_stations.shp")
preprocessed_path = file.path(resultpath, "targets/preprocessed")

<<<<<<< HEAD
if (!file.exists(resultpath)) dir.create(resultpath, recursive = TRUE)
if(!file.exists(preprocessed_path)) dir.create(preprocessed_path, recursive = TRUE)

# define the start and end dates of the period
start_date = as.Date("1981-01-01")
end_date = as.Date("2019-12-31")
=======
resultspath = file.path(rootpath, "results")
check_plot_path = file.path(resultspath, 'ts_plots')
preprocessed_path <- file.path(resultspath, "preprocessed")

lapply(list(resultspath, check_plot_path, preprocessed_path), function(path) {
  if (dir.exists(path)) {
    print("Folder exists!")
  } else {
    dir.create(path)
  }
})

>>>>>>> origin/dev_merged
# ----------------------- pre-processing: select interested European stations -----------------
plan_preprocess = tar_plan(
  tar_target(
    name = 'spanish_rawfiles',
    download_spanish_stations(file.path(rootpath, "Data", "Spain_SAIH"))
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
<<<<<<< HEAD
        name = "write_files",
        command = {
          data.table::fwrite(
            corsica_stations, file = file.path(resultpath,"targets/preprocessed/corsica_stations.csv"))
          data.table::fwrite(
            grdc_stations$output_dt, file = file.path(resultpath,"targets/preprocessed/grdc_updated_stations.csv"))
          data.table::fwrite(
            grdc_old_stations, file = file.path(resultpath, "targets/preprocessed/grdc_old_139stations.csv"))
          data.table::fwrite(
            gsim_stations, file = file.path(resultpath, "targets/preprocessed/gsim_stations.csv"))
          data.table::fwrite(
            rbis_stations, file = file.path(resultpath, "targets/preprocessed/rbis_stations.csv"))
          data.table::fwrite(
            smires_stations$output_ts, file = file.path(resultpath, "targets/preprocessed/smires_stations.csv"))
          data.table::fwrite(
            emr_stations, file = file.path(resultpath, "targets/preprocessed/emr_stations.csv"))
          data.table::fwrite(
            ispra_stations, file = file.path(resultpath,"targets/preprocessed/ispra_stations.csv"))
          data.table::fwrite(
            arpal_stations, file = file.path(resultpath, "targets/preprocessed/arpal_stations.csv"))
          data.table::fwrite(
            arpas_stations, file = file.path(resultpath, "targets/preprocessed/arpas_stations.csv"))
        }, format = "file"
=======
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
  
  tar_target(
    ts_plots,
    lapply(list(corsica_stations$output_ts,
                grdc_stations$output_dt,
                grdc_old_stations,
                smires_stations$output_ts,
                emr_stations$output_ts,
                ispra_stations$output_ts,
                arpal_stations$output_ts,
                arpas_stations$output_ts),
    function(in_dt) {
      lapply(unique(in_dt$gaugeid), 
             function(id) plot_daily_q_timeseries(in_dt[gaugeid == id,], 
                                                  outdir=check_plot_path))
    }
    )
>>>>>>> origin/dev_merged
  )
)
# ----------------------- pre-processing: computing the predictors -----------------
plan_compute_predictors = tar_plan(
  tar_target(
    name = "hr_path",
    command = file.path(datapath,
                        "predictors/HR/waterGAP_streamflow_stations.csv"
    ), format = "file"
  ),
  tar_target(
    name = "hr_predictors",
    command = calc_highres_predictors(hr_path)
  )
)

