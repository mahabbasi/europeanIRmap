grdc_txt2csv <- function(path_list, 
                         start_date ="1981-01-01",
                         end_date ="2019-12-31"){
  
  all_dt <- data.table(dates = seq.Date(as.Date("1901-01-01"),
                                        as.Date("2022-12-31"),
                                        "day"))
  out <- lapply(seq_along(path_list), function(i){
    #extract GRDC unique ID by formatting path
    gaugeno <- strsplit(basename(path_list[i]), '[.]')[[1]][1]
    
    gaugeid <- gaugeno %>% strsplit(., "_") %>% 
      .[[1]] %>% .[1]
    
    # read the GRDC text file for each station
    m <- fread(path_list[i], header = T, sep=";",
               colClasses = c('Date', 'character', 'numeric')) %>%
      setnames('YYYY-MM-DD', 'dates') %>% 
      setnames("Value", gaugeid) %>% 
      .[, - "hh:mm"]
    if (i %% 100 == 0) {
      print (i)
    }
    
    # join the values of station by matching dates
    all_dt <- dplyr::left_join(all_dt, m, by = "dates")
    return(all_dt)
  }) %>% 
    do.call("cbind",.) %>% 
    .[,-"dates"]
  
  output_dt <- cbind(all_dt, out) %>%
    dplyr::na_if(., -999) %>%
    filter(dates >= as.Date(start_date) & dates <= as.Date(end_date))
  
  gaugeids <- output_dt %>% colnames(.) %>% .[-1]
  # find the stations less than 36 months of values
  station_less36_mon <- output_dt %>% 
    dplyr::group_by(month = format(dates, "%Y-%m")) %>% 
    dplyr::summarise_each(funs(mean(.)), gaugeids) %>% 
    tidyr::gather(id, value, gaugeids) %>% 
    group_by(id) %>% 
    dplyr::tally(value >= 0) %>% 
    filter(n < 37)
  # remove the stations less than 36 record
  output_dt_final <- output_dt %>% 
    dplyr::select(-station_less36_mon$id)
  
  gauge_ids <- output_dt_final %>%
    colnames(.) %>% 
    .[-1]
  
  out <- list(
    output_ts = output_dt_final,
    gauge_ids = gauge_ids
  )
  
  return(out)
}

#' select the grdc stations within Europe
#' 
#' @description this function reads the GRDC metadata and excludes the stations outside of
#' the given boundary. the output of the function will be a simple feature `sf` file format.
#' 
grdc_csv2shp <- function(data_path, boundary_path, gauge_ids){
  
  grdc_dt <- data.table::fread(data_path)
  eu_countries <- sf::st_read(dsn = boundary_path)
  joined_sf <- grdc_dt %>% 
    filter(grdc_no %in% gauge_ids) %>%
    sf::st_as_sf(.,
                 coords = c("long", "lat")) %>% 
    st_set_crs(st_crs(eu_countries)) %>% 
    sf::st_intersection(., eu_countries)
  
  lon_lat <- joined_sf %>% 
    sf::st_coordinates() %>% 
    as.data.table()
  
  eu_final_sf <- joined_sf %>% 
    sf::st_drop_geometry() %>% 
    dplyr::select(-c(23:79)) %>% 
    dplyr::mutate(lon_lat) %>% 
    sf::st_as_sf(.,
                 coords = c("X", "Y")) %>% 
    st_set_crs(st_crs(eu_countries))
  
}

select_dataset = function(shp_path, dataset_name = NULL, col_names = NULL){
  m = sf::st_read(dsn = shp_path) %>% 
    dplyr::filter(gaugdts == dataset_name) %>% 
    dplyr::pull(col_names)
}

select_new_grdc_stations = function(path, shp_path){
  # convert and combine all of the updated GRDC stations into a CSV file 
  grdc_updated_stations_dt = grdc_txt2csv(path_list = path)
  
  # find the GRDC stations within the shapefile -> 2110 out of 3484
  interested_grdc_stations = select_dataset(shp_path = shp_path,
                                            dataset_name = "GRDC",
                                            col_names = "gauge_d")
  # get the grdc station id within updated and old datasets.
  updated_grdc_stations_id = intersect(interested_grdc_stations, 
                                       grdc_updated_stations_dt$gauge_ids)
  old_grdc_stations_id = base::setdiff(interested_grdc_stations,
                                       grdc_updated_stations_dt$gauge_ids)
  
  output_dt = grdc_updated_stations_dt$output_ts %>% 
    dplyr::select(one_of(c("dates", updated_grdc_stations_id)))
  
  return(list(old_grdc_stations_id = old_grdc_stations_id,
              updated_grdc_stations_id = updated_grdc_stations_id,
              output_dt = output_dt))
  
}

select_old_grdc_stations = function(path, shp_path, old_grdc_stations_id){
  grdc_dd_id_old = sf::st_read(shp_path) %>% 
    dplyr::filter(gauge_d %in% old_grdc_stations_id) %>% 
    dplyr::pull(dd_id)
  
  output = data.table::fread(path) %>% 
    dplyr::select(one_of(c("date", grdc_dd_id_old))) %>%
    dplyr::mutate(date = as.Date(date, format = "%m/%d/%Y")) %>%
    dplyr::filter(date >= as.Date("1981-01-01")) %>% 
    rename(dates = date)
  
  return(output)
}

select_gsim_stations = function(path, shp_path, dataset_name = "GSIM"){
  
  # find the GSIM stations within the shapefile
  gsim_stations_id = select_dataset(shp_path = shp_path,
                                    dataset_name = dataset_name,
                                    col_names = "dd_id")
  
  output = data.table::fread(path) %>% 
    dplyr::select(one_of(c("date", gsim_stations_id))) %>%
    dplyr::mutate(date = as.Date(date, format = "%m/%d/%Y")) %>%
    dplyr::filter(date >= as.Date("1981-01-01")) %>% 
    rename(dates = date)
  
  return(output)
}

select_smires_stations = function(path, shp_path,
                                  start_date ="1981-01-01",
                                  end_date ="2019-12-31"){
  
  all_dt <- data.table(dates = seq.Date(as.Date("1901-01-01"),
                                        as.Date("2022-12-31"),
                                        "day"))
  smires_stations_id = select_dataset(shp_path = shp_path,
                                      dataset_name = "smires",
                                      col_names = "dd_id")
  out <- lapply(seq_along(path), function(i){
    #extract GRDC unique ID by formatting path
    gaugeno <- strsplit(basename(path[i]), '[.]')[[1]][1]
    
    # read the GRDC text file for each station
    m <- data.table::fread(path[i], header = T, sep=",") %>%
      setnames('time', 'dates') %>% 
      setnames("streamflow", gaugeno) %>% 
      mutate(dates = as.Date(dates))
    
    # join the values of station by matching dates
    all_dt <- dplyr::left_join(all_dt, m, by = "dates")
    return(all_dt)
  }) %>% 
    do.call("cbind",.) %>% 
    .[,-"dates"]
  
  output_dt <- cbind(all_dt, out) %>%
    filter(dates >= as.Date(start_date) & dates <= as.Date(end_date)) %>% 
    dplyr::select(one_of(c("dates", smires_stations_id)))
  
  gaugeids <- output_dt %>% colnames(.) %>% .[-1]
  # find the stations less than 36 months of values
  station_less36_mon <- output_dt %>% 
    dplyr::group_by(month = format(dates, "%Y-%m")) %>% 
    dplyr::summarise_each(funs(mean(.)), gaugeids) %>% 
    tidyr::gather(id, value, gaugeids) %>% 
    group_by(id) %>% 
    dplyr::tally(value >= 0) %>% 
    filter(n < 36)
  # remove the stations less than 36 record
  output_dt_final <- output_dt %>% 
    dplyr::select(-station_less36_mon$id)
  
  gauge_ids <- output_dt_final %>%
    colnames(.) %>% 
    .[-1]
  
  out <- list(
    output_ts = output_dt_final,
    gauge_ids = gauge_ids
  )
  
  return(out)
}
select_corsica_stations = function(path, shp_path){
  
  all_dt <- data.table(dates = seq.Date(as.Date("1901-01-01"),
                                        as.Date("2023-04-30"),
                                        "day"))
  
  # find the corsica stations within the shapefile
  corsica_stations_id = select_dataset(shp_path = shp_path,
                                       dataset_name = "Corsica",
                                       col_names = "gauge_d")
  
  out <- lapply(seq_along(path), function(i){
    #extract GRDC unique ID by formatting path
    gaugeno <- strsplit(basename(path[i]), '[.]')[[1]][1]
    
    gaugeid <- gaugeno %>% strsplit(., "_") %>% 
      .[[1]] %>% .[2]
    
    # read the GRDC text file for each station
    m <- fread(path[i], header = T) %>%
      setnames('date_obs_elab', 'dates') %>% 
      setnames("resultat_obs_elab", gaugeid) %>% 
      dplyr::select(-c("V1", "code_station", "code_site", "date_prod", "code_statut", "libelle_statut",
                       "code_methode", "libelle_methode", "code_qualification", 
                       "libelle_qualification", "longitude", "latitude")) %>% 
      dplyr::mutate(dates = as.Date(dates))
    # join the values of station by matching dates
    all_dt <- dplyr::left_join(all_dt, m, by = "dates")
    return(all_dt)
  }) %>% 
    do.call("cbind",.) %>% 
    .[,-"dates"]
  
  output_dt <- cbind(all_dt, out) %>%
    dplyr::na_if(., -999) %>%
    filter(dates >= as.Date("1981-01-01") & dates <= as.Date("2019-12-31")) %>% 
    dplyr::select(one_of("dates", corsica_stations_id))
  
  gaugeids <- output_dt %>% colnames(.) %>% .[-1]
  # find the stations less than 36 months of values
  station_less36_mon <- output_dt %>% 
    dplyr::group_by(month = format(dates, "%Y-%m")) %>% 
    dplyr::summarise_each(funs(mean(.)), gaugeids) %>% 
    tidyr::gather(id, value, gaugeids) %>% 
    group_by(id) %>% 
    dplyr::tally(value >= 0) %>% 
    filter(n < 36)
  # remove the stations less than 36 record
  output_dt_final <- output_dt %>% 
    dplyr::select(-station_less36_mon$id) %>% 
    mutate(across(-dates, ~ . / 1000))
  
}

select_italian_emr_stations = function(path, shp_path){
  
  all_dt <- data.table(dates = seq.Date(as.Date("1901-01-01"),
                                        as.Date("2023-04-30"),
                                        "day"))
  
  # find the corsica stations within the shapefile
  emr_stations_id = select_dataset(shp_path = shp_path,
                                   dataset_name = "arpae",
                                   col_names = "name")
  
  emr_stations_id = emr_stations_id[! emr_stations_id %in% "emr:Salsominore"]
  emr_stations_id = c(emr_stations_id, "emr:Salsominiore")
  out <- lapply(seq_along(path), function(i){
    #extract emr unique ID by formatting path
    gaugeno = gsub(".csv", "", basename(path[i]))
    # gaugeno = strsplit(basename(path[i]), '[.]')[[1]][1]
    
    gaugeid = gaugeno %>% gsub("emr_", "emr:", .) %>% 
      gsub("_", " ", .)
    
    m = fread(path[i], header = T, skip = 14) %>% .[,-2] %>% 
      rename(dates = "Inizio validità (UTC)")  %>% 
      setnames("Portata media giornaliera (M**3/S)", gaugeid) %>% 
      dplyr::mutate(dates = as.Date(dates))
    # join the values of station by matching dates
    all_dt <- dplyr::left_join(all_dt, m, by = "dates")
    return(all_dt)
  }) %>% 
    do.call("cbind",.) %>% 
    .[,-"dates"]
  
  output_dt <- cbind(all_dt, out) %>%
    # dplyr::na_if(., -999) %>%
    filter(dates >= as.Date("1981-01-01") & dates <= as.Date("2019-12-31")) %>% 
    dplyr::select(one_of("dates", emr_stations_id))
}

select_italian_ispra_stations = function(path, shp_path){
  
  all_dt <- data.table(dates = seq.Date(as.Date("1901-01-01"),
                                        as.Date("2022-12-31"),
                                        "day"))
  
  # find the corsica stations within the shapefile
  ispra_stations_id = select_dataset(shp_path = shp_path,
                                     dataset_name = "ispra",
                                     col_names = "name")
  
  out <- lapply(seq_along(path), function(i){
    #extract emr unique ID by formatting path
    gaugeno = gsub(".csv", "", basename(path[i]))
    
    gaugeid = gaugeno %>%
      gsub("_", ":", .)
    
    m = fread(path[i], header = T, skip = 15) %>% 
      rename(dates = "flowvalue") %>% 
      setnames("dateTime", gaugeid)
    # join the values of station by matching dates
    all_dt <- dplyr::left_join(all_dt, m, by = "dates")
    return(all_dt)
  }) %>% 
    do.call("cbind",.) %>% 
    .[, -"dates"]
  
  output_dt <- cbind(all_dt, out) %>%
    # dplyr::na_if(., -999) %>%
    filter(dates >= as.Date("1981-01-01") & dates <= as.Date("2019-12-31")) %>% 
    dplyr::select(one_of("dates", ispra_stations_id))
}

select_italian_arpal_stations = function(path, shp_path){
  
  all_dt <- data.table(dates = seq.Date(as.Date("1901-01-01"),
                                        as.Date("2022-12-31"),
                                        "day"))
  
  # find the corsica stations within the shapefile
  arpal_stations_id = select_dataset(shp_path = shp_path,
                                     dataset_name = "arpal",
                                     col_names = "name")
  
  out <- lapply(seq_along(path), function(i){
    #extract emr unique ID by formatting path
    gaugeno = gsub(".csv", "", basename(path[i]))
    
    gaugeid = gaugeno
    
    m = fread(path[i], header = TRUE) %>% .[, -1] %>% 
      mutate(date = as.Date(date, format = "%Y_%m_%d"),
             Q = suppressWarnings(as.numeric(Q))) %>% 
      rename( dates = "date") %>% 
      setnames("Q", gaugeid)
    # join the values of station by matching dates
    all_dt <- dplyr::left_join(all_dt, m, by = "dates")
    return(all_dt)
  }) %>% 
    do.call("cbind",.) %>% 
    .[, -"dates"]
  
  output_dt <- cbind(all_dt, out) %>%
    # dplyr::na_if(., -999) %>%
    filter(dates >= as.Date("1981-01-01") & dates <= as.Date("2019-12-31")) %>% 
    dplyr::select(one_of("dates", arpal_stations_id))
}

select_italian_arpas_stations = function(path, shp_path){
  
  all_dt <- data.table(dates = seq.Date(as.Date("1901-01-01"),
                                        as.Date("2022-12-31"),
                                        "day"))
  
  # find the corsica stations within the shapefile
  arpas_stations_id = select_dataset(shp_path = shp_path,
                                     dataset_name = "arpas",
                                     col_names = "name")
  
  out <- lapply(seq_along(path), function(i){
    #extract emr unique ID by formatting path
    gaugeid = gsub(".csv", "", basename(path[i])) %>%
      gsub("_",".", .) %>% 
      paste0("sar:", .)
    
    m = fread(path[i], header = TRUE, skip = 3) %>% 
      dplyr::select(one_of(c("date", "Q"))) %>% 
      mutate(date = as.Date(date, format = "%Y_%m_%d"),
             Q = suppressWarnings(as.numeric(sub(",", ".", Q)))) %>% 
      rename(dates = "date") %>% 
      setnames("Q", gaugeid)
    
    # join the values of station by matching dates
    all_dt <- dplyr::left_join(all_dt, m, by = "dates")
    
    return(all_dt)
  }) %>% 
    do.call("cbind",.) %>% 
    .[, -"dates"]
  
  output_dt <- cbind(all_dt, out) %>%
    # dplyr::na_if(., -999) %>%
    filter(dates >= as.Date("1981-01-01") & dates <= as.Date("2019-12-31")) %>% 
    dplyr::select(one_of("dates", arpas_stations_id))
}
