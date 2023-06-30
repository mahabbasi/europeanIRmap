#---------------------------- UTILITY FUNCTIONS --------------------------------

grdc_txt2csv <- function(path_list, 
                         start_date ="1981-01-01",
                         end_date ="2019-12-31"){
  
  
  out <- lapply(seq_along(path_list), function(i){
    
    #extract GRDC unique ID by formatting path
    gaugeno <- strsplit(basename(path_list[i]), '[.]')[[1]][1]
    gaugeid <- strsplit(gaugeno, "_")[[1]][1]
    
    # read the GRDC text file for each station
    m <- fread(path_list[i], header = T, sep=";",
               colClasses = c('Date', 'character', 'numeric')) %>%
      setnames('YYYY-MM-DD', 'dates') %>% 
      setnames("Value", 'value') %>% 
      .[, - "hh:mm"]
    if (i %% 100 == 0) {
      print (i)
    }
    
    all_dt <- data.table(dates = seq.Date(
      as.Date(paste0(format(min(m$dates, na.rm=T), '%Y'),'-01-01')),
      as.Date(paste0(format(max(m$dates, na.rm=T), '%Y'),'-12-31')),
      "day"))
    
    # join the values of station by matching dates
    all_dt <- merge(all_dt, m, by="dates", all.x=T) %>%
      .[dates >= as.Date(start_date) & dates <= as.Date(end_date),] %>%
      .[, gaugeid := gaugeid]
    
    return(all_dt)
  }) %>% 
    rbindlist
  
  out[value == -999, value := NA] 
  out[, month := format(dates, "%Y-%m")]
  
  # find the stations with records for at least 36 months
  station_less36_mon <- out[gaugeid %in%
                              out[, length(unique(month)), 
                                  by=gaugeid][V1>=36,][[1]]
                            ,
  ]
  
  out <- list(
    output_ts = station_less36_mon,
    gauge_ids = station_less36_mon[, unique(gaugeid)]
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
    st_set_crs(4236) %>% 
    st_transform(st_crs(eu_countries)) %>%
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

select_dataset <- function(shp_path, dataset_name = NULL, col_names = NULL){
  m <- sf::st_read(dsn = shp_path) %>% 
    dplyr::filter(gaugdts == dataset_name) %>% 
    dplyr::pull(col_names)
}

#---------------------------- WORKFLOW FUNCTIONS -------------------------------
download_spanish_stations <- function(output_dir) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  basins <- c('GALICIA%20COSTA','CANTABRICO', 'DUERO', 'EBRO', 'GUADALQUIVIR',
              'GUADIANA', 'JUCAR', 'MI%C3%91O-SIL', 'SEGURA', 'TAJO')
  
  daily_q_paths <- lapply(basins, function(ba) {
    url <- paste0(
      "https://ceh-flumen64.cedex.es/anuarioaforos//anuario-2019-2020//",
      ba,
      "//afliq.csv")
    
    daily_tab_path <- file.path(output_dir, paste0('afliq_', ba, '.csv'))
    
    if (!file.exists(daily_tab_path)) {
      print(paste('Downloading', url))
      download.file(url, daily_tab_path)
    }
    
    return(daily_tab_path)
  }) 
  
  stations_metadata <- lapply(basins, function(ba) {
    url <- paste0(
      "https://ceh-flumen64.cedex.es/anuarioaforos//anuario-2019-2020//",
      ba,
      "//afliq.csv")
    
    daily_tab_path <- file.path(output_dir, paste0('estaf_', ba, '.csv'))
    
    if (!file.exists(daily_tab_path)) {
      print(paste('Downloading', url))
      download.file(url, daily_tab_path)
    }
    
    out_tab <- fread(daily_tab_path, sep=";")
    return(out_tab)
  }) %>%
    rbindlist
  
  
  
  return(list(paths=daily_q_paths,
              metadata=stations_metadata))
}


select_new_grdc_stations <- function(path, shp_path){
  # convert and combine all of the updated GRDC stations into a CSV file 
  grdc_updated_stations_dt <- grdc_txt2csv(path_list = path)
  
  # find the GRDC stations within the shapefile -> 2110 out of 3484
  interested_grdc_stations <- select_dataset(shp_path = shp_path,
                                             dataset_name = "GRDC",
                                             col_names = "gauge_d")
  # get the grdc station id within updated and old datasets.
  output_dt <- grdc_updated_stations_dt$output_ts[
    gaugeid %in% interested_grdc_stations,]
  
  updated_grdc_stations_id <- unique(output_dt$gaugeids)
  
  old_grdc_stations_id <- base::setdiff(interested_grdc_stations,
                                        updated_grdc_stations_id)
  
  return(list(old_grdc_stations_id = old_grdc_stations_id,
              updated_grdc_stations_id = updated_grdc_stations_id,
              output_dt = output_dt))
}

select_old_grdc_stations <- function(path, shp_path, old_grdc_stations_id){
  grdc_dd_id_old <- sf::st_read(shp_path) %>% 
    dplyr::filter(gauge_d %in% old_grdc_stations_id) %>% 
    dplyr::pull(dd_id)
  
  output <- data.table::fread(path) %>% 
    dplyr::select(one_of(c("date", grdc_dd_id_old))) %>%
    dplyr::mutate(date = as.Date(date, format = "%m/%d/%Y")) %>%
    dplyr::filter(date >= as.Date("1981-01-01")) %>% 
    rename(dates = date) %>%
    melt(id.vars='dates', variable.name='gaugeid')
  
  return(output)
}

select_gsim_stations <- function(path, shp_path, dataset_name = "GSIM"){
  
  # find the GSIM stations within the shapefile
  gsim_stations_id <- select_dataset(shp_path = shp_path,
                                     dataset_name = dataset_name,
                                     col_names = "dd_id")
  
  output <- data.table::fread(path) %>% 
    dplyr::select(one_of(c("date", gsim_stations_id))) %>%
    dplyr::mutate(date = as.Date(date, format = "%m/%d/%Y")) %>%
    dplyr::filter(date >= as.Date("1981-01-01")) %>% 
    rename(dates = date)
  
  return(output)
}

select_smires_stations <- function(path, shp_path,
                                   start_date ="1981-01-01",
                                   end_date ="2019-12-31"){
  
  smires_stations_id <- select_dataset(shp_path = shp_path,
                                       dataset_name = "smires",
                                       col_names = "dd_id")
  out <- lapply(seq_along(path), function(i){
    #extract GRDC unique ID by formatting path
    gaugeno <- strsplit(basename(path[i]), '[.]')[[1]][1]
    gaugeid <- gaugeno
    
    
    # read the GRDC text file for each station
    m <- data.table::fread(path[i], header = T, sep=",") %>%
      setnames('time', 'dates') %>% 
      setnames("streamflow", 'value') %>% 
      mutate(dates = as.Date(dates))
    
    # join the values of station by matching dates
    all_dt <- data.table(dates = seq.Date(
      as.Date(paste0(format(min(m$dates, na.rm=T), '%Y'),'-01-01')),
      as.Date(paste0(format(max(m$dates, na.rm=T), '%Y'),'-12-31')),
      "day"))
    
    # join the values of station by matching dates
    all_dt <- merge(all_dt, m, by="dates", all.x=T) %>%
      .[dates >= as.Date(start_date) & dates <= as.Date(end_date),] %>%
      .[, gaugeid := gaugeid]
    
    return(all_dt)
  }) %>% 
    rbindlist
  
  out[value == -999, value := NA] 
  out[, month := format(dates, "%Y-%m")]
  
  # find the stations with records for at least 36 months
  station_less36_mon <- out[gaugeid %in%
                              out[, length(unique(month)), 
                                  by=gaugeid][V1>=36,][[1]]
                            ,
  ]
  
  out <- list(
    output_ts = station_less36_mon,
    gauge_ids = station_less36_mon[, unique(gaugeid)]
  )
  
  return(out)
}

select_corsica_stations <- function(path, shp_path,
                                    start_date ="1981-01-01",
                                    end_date ="2019-12-31"){
  
  # find the corsica stations within the shapefile
  corsica_stations_id <- select_dataset(shp_path = shp_path,
                                        dataset_name = "Corsica",
                                        col_names = "gauge_d")
  
  out <- lapply(seq_along(path), function(i){
    #extract unique ID by formatting path
    gaugeno <- strsplit(basename(path[i]), '[.]')[[1]][1]
    
    gaugeid <- strsplit(gaugeno, "_")[[1]][2]
    
    # read the text file for each station
    m <- fread(path[i], header = T) %>%
      setnames('date_obs_elab', 'dates') %>% 
      setnames("resultat_obs_elab", 'value') %>% 
      dplyr::select(-c("V1", "code_station", "code_site", "date_prod", "code_statut", "libelle_statut",
                       "code_methode", "libelle_methode", "code_qualification", 
                       "libelle_qualification", "longitude", "latitude")) %>% 
      dplyr::mutate(dates = as.Date(dates))
    
    all_dt <- data.table(dates = seq.Date(
      as.Date(paste0(format(min(m$dates, na.rm=T), '%Y'),'-01-01')),
      as.Date(paste0(format(max(m$dates, na.rm=T), '%Y'),'-12-31')),
      "day"))
    
    # join the values of station by matching dates
    all_dt <- merge(all_dt, m, by="dates", all.x=T) %>%
      .[dates >= as.Date(start_date) & dates <= as.Date(end_date),] %>%
      .[, gaugeid := gaugeid]
    
    return(all_dt)
  })  %>% 
    rbindlist
  
  out[value == -999, value := NA] 
  out[, value := value/1000]
  out[, month := format(dates, "%Y-%m")]
  
  # find the stations with records for at least 36 months
  station_less36_mon <- out[gaugeid %in%
                              out[, length(unique(month)), 
                                  by=gaugeid][V1>=36,][[1]]
                            ,
  ]
  
  out <- list(
    output_ts = station_less36_mon,
    gauge_ids = station_less36_mon[, unique(gaugeid)]
  )
  
  return(out)
}

select_italian_emr_stations <- function(path, shp_path,
                                        start_date ="1981-01-01",
                                        end_date ="2019-12-31"){
  
  # find the corsica stations within the shapefile
  emr_stations_id <- select_dataset(shp_path = shp_path,
                                    dataset_name = "arpae",
                                    col_names = "name")
  
  emr_stations_id <- emr_stations_id[! emr_stations_id %in% "emr:Salsominore"]
  emr_stations_id <- c(emr_stations_id, "emr:Salsominiore")
  out <- lapply(seq_along(path), function(i){
    #extract emr unique ID by formatting path
    gaugeno <- gsub(".csv", "", basename(path[i]))
    # gaugeno <- strsplit(basename(path[i]), '[.]')[[1]][1]
    
    gaugeid <- gaugeno %>% gsub("emr_", "emr:", .) %>% 
      gsub("_", " ", .)
    
    m <- fread(path[i], header = T, skip = 14) %>% .[,-2] %>% 
      rename(dates = "Inizio validità (UTC)")  %>% 
      setnames("Portata media giornaliera (M**3/S)", 'value') %>% 
      dplyr::mutate(dates = as.Date(dates))
    
    # join the values of station by matching dates
    all_dt <- data.table(dates = seq.Date(
      as.Date(paste0(format(min(m$dates, na.rm=T), '%Y'),'-01-01')),
      as.Date(paste0(format(max(m$dates, na.rm=T), '%Y'),'-12-31')),
      "day"))
    
    # join the values of station by matching dates
    all_dt <- merge(all_dt, m, by="dates", all.x=T) %>%
      .[dates >= as.Date(start_date) & dates <= as.Date(end_date),] %>%
      .[, gaugeid := gaugeid]
    
    return(all_dt)
  }) %>% 
    rbindlist
  
  out[, month := format(dates, "%Y-%m")]
  
  # find the stations with records for at least 36 months
  station_less36_mon <- out[gaugeid %in%
                              out[, length(unique(month)), 
                                  by=gaugeid][V1>=36,][[1]]
                            ,
  ]
  
  out <- list(
    output_ts = station_less36_mon,
    gauge_ids = station_less36_mon[, unique(gaugeid)]
  )
  
  return(out)
}

select_italian_ispra_stations <- function(path, shp_path,
                                          start_date ="1981-01-01",
                                          end_date ="2019-12-31"){
  
  # find the corsica stations within the shapefile
  ispra_stations_id <- select_dataset(shp_path = shp_path,
                                      dataset_name = "ispra",
                                      col_names = "name")
  
  out <- lapply(seq_along(path), function(i){
    #extract emr unique ID by formatting path
    gaugeno <- gsub(".csv", "", basename(path[i]))
    
    gaugeid <- gaugeno %>%
      gsub("_", ":", .)
    
    m <- fread(path[i], header = T, skip = 15) %>% 
      rename(dates = "flowvalue") %>% 
      setnames("dateTime", 'value')
    
    # join the values of station by matching dates
    all_dt <- data.table(dates = seq.Date(
      as.Date(paste0(format(min(m$dates, na.rm=T), '%Y'),'-01-01')),
      as.Date(paste0(format(max(m$dates, na.rm=T), '%Y'),'-12-31')),
      "day"))
    
    # join the values of station by matching dates
    all_dt <- merge(all_dt, m, by="dates", all.x=T) %>%
      .[dates >= as.Date(start_date) & dates <= as.Date(end_date),] %>%
      .[, gaugeid := gaugeid]
    
    return(all_dt)
  }) %>% 
    rbindlist
  
  out[, month := format(dates, "%Y-%m")]
  
  # find the stations with records for at least 36 months
  station_less36_mon <- out[gaugeid %in%
                              out[, length(unique(month)), 
                                  by=gaugeid][V1>=36,][[1]]
                            ,
  ]
  
  out <- list(
    output_ts = station_less36_mon,
    gauge_ids = station_less36_mon[, unique(gaugeid)]
  )
  
  return(out)
}

select_italian_arpal_stations <- function(path, shp_path,
                                          start_date ="1981-01-01",
                                          end_date ="2019-12-31"){
  
  # find the corsica stations within the shapefile
  arpal_stations_id <- select_dataset(shp_path = shp_path,
                                      dataset_name = "arpal",
                                      col_names = "name")
  
  out <- lapply(seq_along(path), function(i){
    #extract emr unique ID by formatting path
    gaugeno <- gsub(".csv", "", basename(path[i]))
    
    gaugeid <- gaugeno
    
    m <- fread(path[i], header = TRUE) %>% .[, -1] %>% 
      mutate(date = as.Date(date, format = "%Y_%m_%d"),
             Q = suppressWarnings(as.numeric(Q))) %>% 
      rename( dates = "date") %>% 
      setnames("Q", 'value')
    
    # join the values of station by matching dates
    all_dt <- data.table(dates = seq.Date(
      as.Date(paste0(format(min(m$dates, na.rm=T), '%Y'),'-01-01')),
      as.Date(paste0(format(max(m$dates, na.rm=T), '%Y'),'-12-31')),
      "day"))
    
    # join the values of station by matching dates
    all_dt <- merge(all_dt, m, by="dates", all.x=T) %>%
      .[dates >= as.Date(start_date) & dates <= as.Date(end_date),] %>%
      .[, gaugeid := gaugeid]
    
    return(all_dt)
  }) %>% 
    rbindlist
  
  out[, month := format(dates, "%Y-%m")]
  
  # find the stations with records for at least 36 months
  station_less36_mon <- out[gaugeid %in%
                              out[, length(unique(month)), 
                                  by=gaugeid][V1>=36,][[1]]
                            ,
  ]
  
  out <- list(
    output_ts = station_less36_mon,
    gauge_ids = station_less36_mon[, unique(gaugeid)]
  )
  
  return(out)
}

select_italian_arpas_stations <- function(path, shp_path,
                                          start_date ="1981-01-01",
                                          end_date ="2019-12-31"){
  
  all_dt <- data.table(dates = seq.Date(as.Date("1901-01-01"),
                                        as.Date("2022-12-31"),
                                        "day"))
  
  # find the corsica stations within the shapefile
  arpas_stations_id <- select_dataset(shp_path = shp_path,
                                      dataset_name = "arpas",
                                      col_names = "name")
  
  out <- lapply(seq_along(path), function(i){
    #extract emr unique ID by formatting path
    gaugeid <- gsub(".csv", "", basename(path[i])) %>%
      gsub("_",".", .) %>% 
      paste0("sar:", .)
    
    m <- fread(path[i], header = TRUE, skip = 3) %>% 
      dplyr::select(one_of(c("date", "Q"))) %>% 
      mutate(date = as.Date(date, format = "%Y_%m_%d"),
             Q = suppressWarnings(as.numeric(sub(",", ".", Q)))) %>% 
      rename(dates = "date") %>% 
      setnames("Q", 'value')
    
    # join the values of station by matching dates
    all_dt <- data.table(dates = seq.Date(
      as.Date(paste0(format(min(m$dates, na.rm=T), '%Y'),'-01-01')),
      as.Date(paste0(format(max(m$dates, na.rm=T), '%Y'),'-12-31')),
      "day"))
    
    # join the values of station by matching dates
    all_dt <- merge(all_dt, m, by="dates", all.x=T) %>%
      .[dates >= as.Date(start_date) & dates <= as.Date(end_date),] %>%
      .[, gaugeid := gaugeid]
    
    return(all_dt)
  }) %>% 
    rbindlist
  
  out[, month := format(dates, "%Y-%m")]
  
  # find the stations with records for at least 36 months
  station_less36_mon <- out[gaugeid %in%
                              out[, length(unique(month)), 
                                  by=gaugeid][V1>=36,][[1]]
                            ,
  ]
  
  out <- list(
    output_ts = station_less36_mon,
    gauge_ids = station_less36_mon[, unique(gaugeid)]
  )
  
  return(out)
}
