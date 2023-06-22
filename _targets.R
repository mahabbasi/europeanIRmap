library(targets)
library(tarchetypes)
tar_option_set(packages = c("sf", "data.table", "tidyverse", "here"),
               memory = "transient", garbage_collection = TRUE, format = "qs")
tar_source("src/pre_processing_functions.R")
# tar_source("src/used_libraries.R")

# define the data path 
rootpath = here::here()
datapath = file.path(rootpath, "Data/targets")
shp_path = file.path(rootpath, "Data/shp/eu_stations/spatial_joint_eu_stations.shp")

preprocessed_path <- file.path(datapath, "preprocessed")

if (dir.exists(preprocessed_path)) {
  print("Folder exists!")
} else {
  dir.create(preprocessed_path)
}
# ----------------------- pre-processing: select interested European stations -----------------
plan_preprocess = tar_plan(
  tar_target(
    name = "input_files",
    command = list(
      list.files(
        path = file.path(datapath, "raw/2023-04-20_08-00_Europe1"),
        pattern = "txt",
        full.names = TRUE
      ),
      list.files(
        path = file.path(datapath, "raw/corsica"),
        pattern = "csv",
        full.names = TRUE
      ),
      file.path(datapath,
                "raw/harmonized_dataset_cleaned_new.csv"
                ),
      list.files(
        path = file.path(datapath, "raw/Smires/uuid"),
        pattern = "csv",
        full.names = TRUE
        ),
      list.files(
        path = file.path(datapath, "raw/Data_IT/emr"),
        pattern = "csv",
        full.names = TRUE
        ),
      list.files(
        path = file.path(datapath, "raw/Data_IT/ispra"),
        pattern = "csv",
        full.names = TRUE
        ),
      list.files(
        path = file.path(datapath, "raw/Data_IT/lom"),
        pattern = "csv",
        full.names = TRUE
        ),
      list.files(
        path = file.path(datapath, "raw/Data_IT/sar"),
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
        name = "write_files",
        command = {
          data.table::fwrite(
            corsica_stations, file = file.path(datapath,"preprocessed/corsica_stations.csv"))
          data.table::fwrite(
            grdc_stations$output_dt, file = file.path(datapath,"preprocessed/grdc_updated_stations.csv"))
          data.table::fwrite(
            grdc_old_stations, file = file.path(datapath, "preprocessed/grdc_old_139stations.csv"))
          data.table::fwrite(
            gsim_stations, file = file.path(datapath, "preprocessed/gsim_stations.csv"))
          data.table::fwrite(
            rbis_stations, file = file.path(datapath, "preprocessed/rbis_stations.csv"))
          data.table::fwrite(
            smires_stations$output_ts, file = file.path(datapath, "preprocessed/smires_stations.csv"))
          data.table::fwrite(
            emr_stations, file = file.path(datapath, "preprocessed/emr_stations.csv"))
          data.table::fwrite(
            ispra_stations, file = file.path(datapath,"preprocessed/ispra_stations.csv"))
          data.table::fwrite(
            arpal_stations, file = file.path(datapath, "preprocessed/arpal_stations.csv"))
          data.table::fwrite(
            arpas_stations, file = file.path(datapath, "preprocessed/arpas_stations.csv"))
        }, format = "file"
  )
)
