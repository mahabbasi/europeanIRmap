
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
      setnames(c('YYYY-MM-DD',"Value"),
               c('dates', 'value')
               ) %>% 
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
  
  out[value %in% c(-999, -99, -9999, 999, 9999), value := NA] 
  
  station_less36_mon <- select_min_fullmonths(out, min_months = 36)

  function_output <- list(
    output_ts = station_less36_mon,
    gauge_ids = station_less36_mon[, unique(gaugeid)]
  )
  
  return(function_output)
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
<<<<<<< HEAD
    st_set_crs(st_crs(4326)) %>% 
=======
    st_set_crs(4236) %>% 
    st_transform(st_crs(eu_countries)) %>%
>>>>>>> origin/dev_merged
    sf::st_intersection(., eu_countries)
  
  lon_lat <- joined_sf %>% 
    sf::st_coordinates() %>% 
    data.table::as.data.table()
  
  eu_final_sf <- joined_sf %>% 
    sf::st_drop_geometry() %>% 
    dplyr::select(-c(23:79)) %>% 
    dplyr::mutate(lon_lat) %>% 
    sf::st_as_sf(.,
                 coords = c("X", "Y")) %>% 
    st_set_crs(st_crs(4326))
  
  return(eu_final_sf)
}

select_dataset <- function(shp_path, dataset_name = NULL, col_names = NULL){
  m <- sf::st_read(dsn = shp_path) %>% 
    dplyr::filter(gaugdts == dataset_name) %>% 
    dplyr::pull(col_names)
}

select_min_fullmonths <- function(in_dt, min_months = 36) {
  #out[value %in% c(-999, -99, -9999, 999, 9999), value := NA] 
  in_dt[, month := format(dates, "%Y-%m")]
  
  in_dt[, missing_days_month := sum(is.na(value)), by=c('gaugeid', 'month')]
  
  # find the stations with records for at least 36 months
  station_less36_mon <- in_dt[gaugeid %in%
                                in_dt[missing_days_month == 0, 
                                      length(unique(month)), 
                                      by=gaugeid][V1>=36,][[1]]
                              ,][missing_days_month == 0,]
  
  return(station_less36_mon[, c('gaugeid', 'dates', 'value'), with=F])
}

#------ flag_outliers ------
#' Flag outliers
#'
#' Flag potential outliers in daily discharge records for a given gauging
#' station following the criteria developed for GSIM by
#' [Gudmundsson et al. (2018)](https://essd.copernicus.org/articles/10/787/2018/).
#'
#' @param in_gaugetab \link[data.table]{data.table} containing formatted daily
#' discharge record from gauging station.
#'
#' @details Criteria to flag a daily discharge value (Qt) as a potential outlier include:
#' \itemize{
#'   \item Negative values (Qt < 0)
#'   \item At least ten identical consecutive discharge values (for Qt > 0)
#'   \item |log(Qt + 0.01) - mean| are larger than the mean values of log(Q + 0.01)
#'   plus or minus 6 times the standard deviation of log(Q + 0.01) computed for
#'   that calendar day for the entire length of the series. The mean and SD are
#'   computed for a 5-day window centred on the calendar day to ensure that a
#'   sufficient amount of data is considered. The log-transformation is used to
#'   account for the skewness of the distribution of daily streamflow values.
#'   \item Qt for which value != Calculated discharge in the record
#' }
#'
#' @return \link[data.table]{data.table} of daily discharge records with additional
#' columns for outlier flags
#'
#' @source Gudmundsson, L., Do, H. X., Leonard, M., & Westra, S. (2018). The Global
#'   Streamflow Indices and Metadata Archive (GSIM) – Part 2: Quality control,
#'   time-series indices and homogeneity assessment. Earth System Science Data,
#'   10(2), 787–804. https://doi.org/10.5194/essd-10-787-2018
#'
#' @export
flag_outliers <- function(in_gaugetab) {
  in_gaugetab[value %in% c(-999, -99, -9999, 999, 9999), value := NA]
  
  in_gaugetab[, `:=`(jday = format(as.Date(dates), '%j'),#Julian day
                     q_rleid = rleid(value),#Identify each group of consecutive values
                     flag_dryver = 0)] #Create flag field)
  
  #Flag negative values
  in_gaugetab[value < 0, flag_dryver := flag_dryver + 1]
  
  #Flag when more than 10 identical values in a row or when a single zero-flow
  in_gaugetab[, flag_dryver := flag_dryver +
                ((value > 0) & (.N > 10)) +
                ((value == 0) & (.N == 1)),
              by=q_rleid]
  
  #Flag |log(Q + 0.01) - mean| > 6SD for julian day mean and SD of 5d mean of log(Q + 0.01)
  in_gaugetab[, logmean5d := frollapply(log(value + 0.01), n = 5, align='center',
                                        FUN=mean, na.rm = T)] %>% #Compute 5-day mean of log(Q+0.01)
    .[, `:=`(jdaymean = mean(logmean5d, na.rm = T),
             jdaysd = sd(logmean5d, na.rm = T)),
      by = jday] %>% #Compute mean and SD of 5-day mean of log(Q + 0.01) by Julian day
    .[abs(log(value + 0.01) - jdaymean) > (6 * jdaysd),
      flag_dryver := flag_dryver + 1]
  
  return(in_gaugetab)
}


#------ plot_daily_q_timeseries ----------------------
#' Plot a daily time series
#'
#' Creates a plot of daily discharge ({m^3}/s) for a gauging station,
#' with flags for 0-flow values and potential outliers.
#' Save plot to png if path is provided.
#'
#' @param gaugestats_record \link[data.table]{data.table} of formatted daily
#' discharge records for a single station. Must contain at least five columns:
#' \code{gaugeid, dates, value, flag_dryver} \cr
#' @param outpath (character) path for writing output png plot (default is no plotting).
#'
#' @details the output graphs show the time series of daily streamflow values
#' for the station. For the flagging criteria, see documentation for \code{\link{flag_outliers}}.
#' \itemize{
#'   \item The y-axis is square-root transformed.
#'   \item Individual points show daily discharge values (in {m^3}/s).
#'   \item blue lines link daily values (which may result in unusual patterns due to missing years).
#'   \item red points are zero-flow flow values.
#'   \item green points are non-zero flow daily values statistically flagged as potential outliers .
#'   \item black points are zero-flow values flagged as potential outliers.
#' }
#'
#' @return plot
#'
#' @export
plot_daily_q_timeseries <- function(in_dt, outdir=NULL) {
  #Read and format discharge records
  if (in_dt[,.N>1] &
      ('value' %in% names(in_dt))) {
    gaugetab <- flag_outliers(copy(in_dt))
  } 
  
  # else {
  #   gaugetab <- readformatGRDC(in_dt$path) %>%
  #     flagGRDCoutliers %>%
  #     .[, dates := as.Date(dates)] %>%
  #     .[!is.na(value), missingdays := diny(year)-.N, by= 'year']
  # }
  
  
  #Plot time series
  qtiles <- union(gaugetab[, min(value, na.rm=T)],
                  gaugetab[, quantile(value, probs=seq(0, 1, 0.1), na.rm=T)])
  
  rawplot <- ggplot(gaugetab,
                    aes(x=dates, y=value)) +
    geom_line(color='#045a8d', size=1, alpha=1/5) +
    geom_point(data = gaugetab[flag_dryver == 0 & value > 0,],
               color='#045a8d', size=1, alpha=1/3) +
    geom_point(data = gaugetab[flag_dryver > 0 & value > 0,],
               color='green') +
    geom_point(data = gaugetab[flag_dryver == 0 & value == 0,],
               color='red') +
    geom_point(data = gaugetab[flag_dryver > 0 & value == 0,],
               color='black') +
    scale_y_sqrt(breaks=qtiles, labels=qtiles) +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    labs(y='Discharge (m3/s)',
         title=paste0('GRDC: ', in_dt$gaugeid)) +
    coord_cartesian(expand=0, clip='off')+
    theme_bw() +
    theme(axis.text.x = element_text(angle = 44, hjust=1),
          axis.text.y = element_text())
  
  # if (showmissing) {
  #   rawplot <- rawplot +
  #     geom_point(data=gaugetab[missingdays >= maxgap,], 
  #                color='black', alpha=1/10)
  # }
  
  if (!is.null(outdir)) {
    outpath <-  file.path(outdir, paste0(unique(in_dt$gaugeid), '.png'))
    
    if (!(file.exists(outpath))) {
      print(paste0("Saving ", outpath, "..."))
      ggsave(filename = outpath,
             plot = rawplot, device = 'png',
             width = 10, height = 10, units='in', dpi = 300)
    }
  } else {
    return(outpath)
  }
}

#------ plotGSIMtimeseries ----------------------
#' Plot a GSIM time series
#'
#' Creates a plot of monthly discharge ({m^3}/s) for a GSIM gauging stations
#' Save plot to png if path is provided.
#'
#' @param GSIMgaugestats_record a data containing structure  with
#' a column called "path" towards a standard montly GSIM text file containing.
#' In this project, e.g. the output from \code{\link{comp_GSIMdurfreq}}.
#' @param outpath (character) path for writing output png plot (default is no plotting).
#' @param maxgap (integer) threshold number of missing daily records to consider a calendar year unfit for analysis.
#' @param showmissing (logical) whether to show records in years with number of missing daily records beyond \code{maxgap}.
#'
#' @details Daily streamflow records from GSIM stations are unavailable. Therefore,
#' the graph shows the following:
#' \itemize{
#'   \item The y-axis is square-root transformed.
#'   \item Blue points: mean monthly discharge
#'   \item Light blue background shading: mean ± 2SD monthly discharge
#'   \item Black points: minimum and maximum monthly discharge
#'   \item Red points show minimum monthly discharge values equal to 0
#'   \item Purple points show months for which all daily discharge values are equal to 0.
#' }
#'
#' @return plot
#'
#' @source Gudmundsson, L., Do, H. X., Leonard, M., & Westra, S. (2018). The Global
#'   Streamflow Indices and Metadata Archive (GSIM) – Part 2: Quality control,
#'   time-series indices and homogeneity assessment. Earth System Science Data,
#'   10(2), 787–804. https://doi.org/10.5194/essd-10-787-2018
#'
#' @export
plot_gsim_q_timeseries <- function(in_dt, outpath=NULL) {
  #Read and format discharge records
  gaugetab <- copy(in_dt) %>%
    .[!is.na(MEAN),] %>%
    setorder(date)
  
  #Format for plotting: compute MEAN - 2*SD if > MIN, otherwise MIN
  gaugetab[, `:=`(ribbonlow = max(c(MEAN-2*SD, MIN)),
                  ribbonhigh = min(c(MEAN+2*SD, MAX))
  ), by=date]
  
  #Plot time series
  qtiles <- gaugetab[, quantile(union(MIN, MAX), probs=seq(0, 1, 0.1))]
  
  rawplot <- ggplot(gaugetab[missingdays < maxgap,], aes(x=date, y=MEAN)) +
    geom_line(color='#045a8d', size=1, alpha=1/3) +
    geom_point(data=gaugetab[(missingdays < maxgap) & (MEAN>0),],
               color='#045a8d', size=1, alpha=1/2) +
    geom_point(data=gaugetab[(missingdays < maxgap) & (MEAN==0),],
               color='red', size=1, alpha=1/2) +
    geom_ribbon(aes(ymin = ribbonlow, ymax = ribbonhigh), color='lightblue', alpha=1/4) +
    geom_point(data=gaugetab[(missingdays < maxgap) & (MIN>0),],
               aes(y = MIN), color='black', alpha=1/2) +
    geom_point(data=gaugetab[(missingdays < maxgap) & (MIN==0),],
               aes(y = MIN), color='darkred', alpha=1/2) +
    geom_point(data=gaugetab[(missingdays < maxgap) & (MAX>0),],
               aes(y = MAX), color='black', alpha=1/2) +
    geom_point(data=gaugetab[(missingdays < maxgap) & (MAX==0),],
               aes(y = MAX), color='purple', alpha=1/2) +
    scale_y_sqrt(breaks=qtiles, labels=qtiles) +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    labs(y='Discharge (m3/s)',
         title=paste0('GSIM: ', GSIMgaugestats_record$gsim_no)) +
    coord_cartesian(expand=0, clip='off')+
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust=1),
          axis.text.y = element_text())
  
  if (showmissing) {
    gaugetab_removed <- gaugetab[missingdays > maxgap,]
    rawplot <- rawplot +
      geom_point(data=gaugetab[MEAN>0,],
                 color='#045a8d', size=1, alpha=1/5) +
      geom_point(data=gaugetab[(missingdays < maxgap) & (MEAN==0),],
                 color='red', size=1, alpha=1/5) +
      geom_point(data=gaugetab[(missingdays < maxgap) & (MIN>0),],
                 aes(y = MIN), color='black', alpha=1/5) +
      geom_point(data=gaugetab[(missingdays < maxgap) & (MIN==0),],
                 aes(y = MIN), color='darkred', alpha=1/5) +
      geom_point(data=gaugetab[(missingdays < maxgap) & (MAX>0),],
                 aes(y = MAX), color='black', alpha=1/5) +
      geom_point(data=gaugetab[(missingdays < maxgap) & (MAX==0),],
                 aes(y = MAX), color='purple', alpha=1/5)
  }
  
  if (!is.null(outdir)) {
    outpath <-  file.path(outdir, paste0(unique(in_dt$gaugeid), '.png'))
    
    if (!(file.exists(outpath))) {
      print(paste0("Saving ", outpath, "..."))
      ggsave(filename = outpath, plot = rawplot, device = 'png',
             width = 10, height = 10, units='in', dpi = 300)
    }
  } else {
    return(rawplot)
  }
}

#---------------------------- WORKFLOW FUNCTIONS -------------------------------
download_spanish_stations <- function(output_dir) {
  #https://ceh.cedex.es/anuarioaforos/demarcaciones.asp
  #indroea: Indicativo oficial de la estación de aforos
  #fecha: dd/mm/yyyy
  #altura: water height - m
  #caudal: discharge - m3/s
  #suprest: drainage area in - km2
  
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
  }) %>% unlist
  
  stations_metadata <- lapply(basins, function(ba) {
    url <- paste0(
      "https://ceh-flumen64.cedex.es/anuarioaforos//anuario-2019-2020//",
      ba,
      "//estaf.csv")
    
    metadata_path <- file.path(output_dir, paste0('estaf_', ba, '.csv'))
    
    if (!file.exists(metadata_path)) {
      print(paste('Downloading', url))
      download.file(url, metadata_path)
    }
    
    out_tab <- fread(metadata_path, sep=";")
    return(out_tab)
  }) %>%
    rbindlist
  
  return(list(paths=daily_q_paths,
              metadata=stations_metadata))
}

select_spanish_stations <- function(in_paths, in_metadata,
                                    start_date ="1981-01-01",
                                    end_date ="2019-12-31") {
  
  out <- lapply(seq_along(in_paths), function(i){

    # read the GRDC text file for each station
    m <- data.table::fread(in_paths[i], header = T, sep=";") %>%
      setnames(c('fecha', 'caudal', 'indroea'),
               c('dates', 'value', 'gaugeid')
               ) %>% 
      mutate(dates = as.Date(dates, format = "%d/%m/%Y"))
    
    # join the values of station by matching dates
    all_dt <- m[,
                list(
                  dates = seq.Date(
                    as.Date(paste0(format(min(dates, na.rm=T), '%Y'),'-01-01')),
                    as.Date(paste0(format(max(dates, na.rm=T), '%Y'),'-12-31')),
                    "day")
                ),
                by=gaugeid
    ]
    
    # join the values of station by matching dates
    all_dt <- merge(all_dt, m, by=c("dates", "gaugeid"), all.x=T) %>%
      .[dates >= as.Date(start_date) & dates <= as.Date(end_date),]
    
    return(all_dt)
  }) %>% 
    rbindlist
  
  station_less36_mon <- select_min_fullmonths(out, min_months = 36) 
  
  function_output <- list(
    output_ts = station_less36_mon,
    gauge_ids = station_less36_mon[, unique(gaugeid)]
  )
  
  return(function_output)
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
  
  updated_grdc_stations_id <- unique(output_dt$gaugeid)
  
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
  
<<<<<<< HEAD
  # read the global variable
  start_date <<- start_date
  
  output = data.table::fread(path) %>% 
    dplyr::select(one_of(c("date", grdc_dd_id_old))) %>%
    dplyr::mutate(date = as.Date(date, format = "%m/%d/%Y")) %>%
    dplyr::filter(date >= start_date) %>% 
    rename(dates = date)
=======
  output <- data.table::fread(path) %>% 
    dplyr::select(one_of(c("date", grdc_dd_id_old))) %>%
    dplyr::mutate(date = as.Date(date, format = "%m/%d/%Y")) %>%
    dplyr::filter(date >= as.Date("1981-01-01")) %>% 
    rename(dates = date) %>%
    melt(id.vars='dates', variable.name='gaugeid')
>>>>>>> origin/dev_merged
  
  return(output)
}

select_gsim_stations <- function(path, shp_path, dataset_name = "GSIM"){
  
  # find the GSIM stations within the shapefile
<<<<<<< HEAD
  gsim_stations_id = select_dataset(shp_path = shp_path,
                                    dataset_name = dataset_name,
                                    col_names = "dd_id")
  # read the global variable
  start_date <<- start_date
=======
  gsim_stations_id <- select_dataset(shp_path = shp_path,
                                     dataset_name = dataset_name,
                                     col_names = "dd_id")
>>>>>>> origin/dev_merged
  
  output <- data.table::fread(path) %>% 
    dplyr::select(one_of(c("date", gsim_stations_id))) %>%
    dplyr::mutate(date = as.Date(date, format = "%m/%d/%Y")) %>%
    dplyr::filter(date >= start_date) %>% 
    rename(dates = date)

  return(output)
}

<<<<<<< HEAD
select_smires_stations = function(path, shp_path){
  
  # read the global variables
  start_date <<- start_date
  end_date <<- end_date
=======
select_smires_stations <- function(path, shp_path,
                                   start_date ="1981-01-01",
                                   end_date ="2019-12-31"){
>>>>>>> origin/dev_merged
  
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
  
<<<<<<< HEAD
  output_dt <- cbind(all_dt, out) %>%
    filter(dates >= start_date & dates <= end_date) %>% 
    dplyr::select(one_of(c("dates", smires_stations_id)))
=======
  out[value %in% c(-999, -99, -9999, 999, 9999), value := NA] 
>>>>>>> origin/dev_merged
  
  # find the stations with records for at least 36 months
  station_less36_mon <- select_min_fullmonths(out, min_months = 36)
  
  function_output <- list(
    output_ts = station_less36_mon,
    gauge_ids = station_less36_mon[, unique(gaugeid)]
  )
  
  return(function_output)
}
<<<<<<< HEAD
select_corsica_stations = function(path, shp_path){
  
  # assign start and end date from the global variables
  start_date <<- start_date
  end_date <<- end_date
  
  all_dt <- data.table(dates = seq.Date(as.Date("1901-01-01"),
                                        as.Date("2023-04-30"),
                                        "day"))
=======

select_corsica_stations <- function(path, shp_path,
                                    start_date ="1981-01-01",
                                    end_date ="2019-12-31"){
>>>>>>> origin/dev_merged
  
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
  
<<<<<<< HEAD
  output_dt <- cbind(all_dt, out) %>%
    dplyr::na_if(., -999) %>%
    filter(dates >= start_date & dates <= end_date) %>% 
    dplyr::select(one_of("dates", corsica_stations_id))
=======
  out[value %in% c(-999, -99, -9999, 999, 9999), value := NA] 
  out[, value := value/1000]

  # find the stations with records for at least 36 months
  station_less36_mon <- select_min_fullmonths(out, min_months = 36)
>>>>>>> origin/dev_merged
  
  function_output <- list(
    output_ts = station_less36_mon,
    gauge_ids = station_less36_mon[, unique(gaugeid)]
  )
  
  return(function_output)
}

<<<<<<< HEAD
select_italian_emr_stations = function(path, shp_path){
  
  # assign start and end date from the global variables
  start_date <<- start_date
  end_date <<- end_date
  
  all_dt <- data.table(dates = seq.Date(as.Date("1901-01-01"),
                                        as.Date("2023-04-30"),
                                        "day"))
=======
select_italian_emr_stations <- function(path, shp_path,
                                        start_date ="1981-01-01",
                                        end_date ="2019-12-31"){
>>>>>>> origin/dev_merged
  
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
  
<<<<<<< HEAD
  output_dt <- cbind(all_dt, out) %>%
    # dplyr::na_if(., -999) %>%
    filter(dates >= start_date & dates <= end_date) %>% 
    dplyr::select(one_of("dates", emr_stations_id))
}

select_italian_ispra_stations = function(path, shp_path){
  
  # assign start and end date from the global variables
  start_date <<- start_date
  end_date <<- end_date
  
  all_dt <- data.table(dates = seq.Date(as.Date("1901-01-01"),
                                        as.Date("2022-12-31"),
                                        "day"))
=======
  # find the stations with records for at least 36 months
  station_less36_mon <- select_min_fullmonths(out, min_months = 36)
  
  function_output <- list(
    output_ts = station_less36_mon,
    gauge_ids = station_less36_mon[, unique(gaugeid)]
  )
  
  return(function_output)
}

select_italian_ispra_stations <- function(path, shp_path,
                                          start_date ="1981-01-01",
                                          end_date ="2019-12-31"){
>>>>>>> origin/dev_merged
  
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
  
<<<<<<< HEAD
  output_dt <- cbind(all_dt, out) %>%
    # dplyr::na_if(., -999) %>%
    filter(dates >= start_date & dates <= end_date) %>% 
    dplyr::select(one_of("dates", ispra_stations_id))
}

select_italian_arpal_stations = function(path, shp_path){
  
  # assign start and end date from the global variables
  start_date <<- start_date
  end_date <<- end_date
  
  all_dt <- data.table(dates = seq.Date(as.Date("1901-01-01"),
                                        as.Date("2022-12-31"),
                                        "day"))
=======
  # find the stations with records for at least 36 months
  station_less36_mon <- select_min_fullmonths(out, min_months = 36)
  
  function_output <- list(
    output_ts = station_less36_mon,
    gauge_ids = station_less36_mon[, unique(gaugeid)]
  )
  
  return(function_output)
}

select_italian_arpal_stations <- function(path, shp_path,
                                          start_date ="1981-01-01",
                                          end_date ="2019-12-31"){
>>>>>>> origin/dev_merged
  
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
  
<<<<<<< HEAD
  output_dt <- cbind(all_dt, out) %>%
    # dplyr::na_if(., -999) %>%
    filter(dates >= start_date & dates <= end_date) %>% 
    dplyr::select(one_of("dates", arpal_stations_id))
=======
  # find the stations with records for at least 36 months
  station_less36_mon <- select_min_fullmonths(out, min_months = 36)
  
  function_output <- list(
    output_ts = station_less36_mon,
    gauge_ids = station_less36_mon[, unique(gaugeid)]
  )
  
  return(function_output)
>>>>>>> origin/dev_merged
}

select_italian_arpas_stations <- function(path, shp_path,
                                          start_date ="1981-01-01",
                                          end_date ="2019-12-31"){
  
  # assign start and end date from the global variables
  start_date <<- start_date
  end_date <<- end_date
  
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
  
<<<<<<< HEAD
  output_dt <- cbind(all_dt, out) %>%
    # dplyr::na_if(., -999) %>%
    filter(dates >= start_date & dates <= end_date) %>% 
    dplyr::select(one_of("dates", arpas_stations_id))
}

### ----------------- Functions for calculating the High resolutions predictors ----------------
#' add_date: this function is used to generate a sequence monthly date within our time span.
#' 
#' @param mat (matrix) the time series of WaterGAP streamflow for all the stations.
#'  
add_date <- function(mat){
  # assign start and end date from the global variables
  end_date <<- end_date
  
  mat %>% 
    as.data.frame(.) %>% 
    replace(is.na(.), 0) %>%
    mutate(date = seq.Date(as.Date("1980-01-01"),
                           end_date,
                           "month"))
}

#' calculating the minimum streamflow of the past three months for the gauging stations.
#' 
#'  @param tbl (vector) the gauging station streamflow vector in either dataframe or tibble class. 
min_past_3_month <- function(tbl){
  tbl %>%
    mutate(Qlag1 = lag(.$Q, 1),
           Qlag2 = lag(.$Q, 2),
           Qlag3 = lag(.$Q, 3)) %>%
    rowwise() %>%
    mutate(min_p3m = min(c(Qlag1, Qlag2, Qlag3)))
}

#' calculating the average streamflow of the past three months for the gauging stations.
#' 
#'  @param tbl (vector) the gauging station streamflow vector in either dataframe or tibble class. 
#'  
mean_past_3_month <- function(tbl){
  tbl %>%
    mutate(Qlag1 = lag(.$Q, 1),
           Qlag2 = lag(.$Q, 2),
           Qlag3 = lag(.$Q, 3)) %>%
    rowwise() %>%
    mutate(mean_p3m = mean(c(Qlag1, Qlag2, Qlag3)))
}

#' calculating the minimum streamflow of the past 12 months for the gauging stations.
#' 
#'  @param tbl (vector) the gauging station streamflow vector in either dataframe or tibble class. 
#'  
min_past_12_month <- function(tbl){
  tbl %>%
    mutate(Qlag1 = lag(.$Q, 1), Qlag2 = lag(.$Q, 2), Qlag3 = lag(.$Q, 3),
           Qlag4 = lag(.$Q, 4), Qlag5 = lag(.$Q, 5), Qlag6 = lag(.$Q, 6),
           Qlag7 = lag(.$Q, 7), Qlag8 = lag(.$Q, 8), Qlag9 = lag(.$Q, 9),
           Qlag10 = lag(.$Q, 10), Qlag11 = lag(.$Q, 11), Qlag12 = lag(.$Q, 12)) %>%
    rowwise() %>%
    mutate(min_p12m = min(c(Qlag1, Qlag2, Qlag3,Qlag4, Qlag5, Qlag6,
                            Qlag7, Qlag8, Qlag9,Qlag10, Qlag11, Qlag12)))
}

#' calculating the average streamflow of the past 12 months for the gauging stations.
#' 
#'  @param tbl (vector) the gauging station streamflow vector in either dataframe or tibble class. 
#'  
mean_past_12_month <- function(tbl){
  tbl %>%
    mutate(Qlag1 = lag(.$Q, 1), Qlag2 = lag(.$Q, 2), Qlag3 = lag(.$Q, 3),
           Qlag4 = lag(.$Q, 4), Qlag5 = lag(.$Q, 5), Qlag6 = lag(.$Q, 6),
           Qlag7 = lag(.$Q, 7), Qlag8 = lag(.$Q, 8), Qlag9 = lag(.$Q, 9),
           Qlag10 = lag(.$Q, 10), Qlag11 = lag(.$Q, 11), Qlag12 = lag(.$Q, 12)) %>%
    rowwise() %>%
    mutate(mean_p12m = mean(c(Qlag1, Qlag2, Qlag3,Qlag4, Qlag5, Qlag6,
                              Qlag7, Qlag8, Qlag9,Qlag10, Qlag11, Qlag12)))
}
#' calculating the standard deviation of streamflow for 12 months for the gauging stations.
#' 
#'  @param tbl (vector) the gauging station streamflow vector in either dataframe or tibble class. 
#'  @param start_date (character) the start date of the time period
#'  
sd_mon <- function(tbl){
  
  # assign start and end date from the global variables
  start_date <<- start_date

  tbl %>%
    filter(date >= start_date) %>%
    as.vector(.) %>% .$Q %>% matrix(., ncol = 12, byrow = TRUE) %>%
    as.data.frame(.) %>% `colnames<-`(month.abb) %>%
    summarise_all(., sd)
} 

#' calculating the average of streamflow for 12 months for the gauging stations.
#' 
#'  @param tbl (vector) the gauging station streamflow vector in either dataframe or tibble class. 
#'  @param start_date (character) the start date of the time period
#'
mean_mon <- function(tbl){
  
  # assign start and end date from the global variables
  start_date <<- start_date

  tbl %>%
    filter(date >= start_date) %>%
    as.vector(.) %>% .$Q %>% matrix(., ncol = 12, byrow = TRUE) %>%
    as.data.frame(.) %>% `colnames<-`(month.abb) %>%
    summarise_all(., mean)
}

#' combining all the above functions to calculate the high-resolution predictors for all the gauging stations.
#' 
#'  @param watergap_raw_path (character) the path for streamflow time series at all the gauging stations derived
#'  from `highres_streamflow_ext`. 
#'  
#'  @return a list with the following elements:
#'  \itemize{
#'    \item min_p3m - A list of the minimum streamflow of the past three months for all the gaugins stations.
#'    \item mean_p3m - A list of the average streamflow of the past three months for all the gaugins stations.
#'    \item min_p12m - A list of the minimum streamflow of the past 12 months for all the gaugins stations.
#'    \item mean_p3m - A list of the average streamflow of the past 12 months for all the gaugins stations.
#'    \item sd - A list of the standard deviation of streamflow of 12 calender months for all the gaugins stations.
#'    \item cv - A list of the Coefficient of variation of streamflow of 12 calender months for all the gaugins stations.
#'  }
#'
calc_highres_predictors <- function(watergap_raw_path){
  
  # assign start and end date from the global variables
  start_date <<- start_date
  end_date <<- end_date
  
  # read the monthly downscaled watergap streamflow for stations
  waterGap_streamflow = data.table::fread(watergap_raw_path)
  stations_dd_id = colnames(waterGap_streamflow)
  # defining the empty lists 
  list_min_p3m = list()
  list_mean_p3m = list()
  list_min_p12m = list()
  list_mean_p12m = list()
  sd_list = list()
  cv_list = list()
  
  # looping over stations to compute the high resolution predictors
  for (i in seq_along(stations_dd_id)) {
    WaterGap_past_3_month_min = waterGap_streamflow %>% 
      add_date(.) %>%
      dplyr::select(date, Q = stations_dd_id[i]) %>%
      min_past_3_month(.)
    WaterGap_past_3_month_mean = waterGap_streamflow %>% 
      add_date(.) %>%
      dplyr::select(date, Q = stations_dd_id[i]) %>%
      mean_past_3_month(.)
    WaterGap_past_12_month_min = waterGap_streamflow %>% 
      add_date(.) %>%
      dplyr::select(date, Q = stations_dd_id[i]) %>%
      min_past_12_month(.)
    WaterGap_past_12_month_mean = waterGap_streamflow %>% 
      add_date(.) %>%
      dplyr::select(date, Q = stations_dd_id[i]) %>%
      mean_past_12_month(.)
    df_mean = waterGap_streamflow %>% 
      add_date(.) %>%
      dplyr::select(date, Q = stations_dd_id[i]) %>% 
      mean_mon(.)
    df_sd = waterGap_streamflow %>% 
      add_date(.) %>%
      dplyr::select(date, Q = stations_dd_id[i]) %>% 
      sd_mon(.)
    cv = df_sd / df_mean
    # repeat the interannual predictors (sd, cv) over the period
    # the period is 40 years.
    sd_ts = rep(as.numeric(df_sd), 40) %>%
      as.data.frame(.) %>% `colnames<-`("sd") %>%
      mutate(date = seq.Date(as.Date("1980-01-01"),
                             end_date, "month"),
             .before = "sd")
    
    cv_ts = rep(as.numeric(cv), 40) %>%
      as.data.frame(.) %>% `colnames<-`("cv") %>%
      mutate(date = seq.Date(as.Date("1980-01-01"),
                             end_date, "month"),
             .before = "cv")
    # store the predictors into the lists
    list_min_p3m[[stations_dd_id[i]]] = WaterGap_past_3_month_min
    list_mean_p3m[[stations_dd_id[i]]] = WaterGap_past_3_month_mean
    list_min_p12m[[stations_dd_id[i]]] = WaterGap_past_12_month_min
    list_mean_p12m[[stations_dd_id[i]]] = WaterGap_past_12_month_mean
    sd_list[[stations_dd_id[i]]] = sd_ts
    cv_list[[stations_dd_id[i]]] = cv_ts
    print(i)
  }
  # convert the lists to df
  cat("the lists of HR predictors are been converting to dataframe.\n")
  # minimum of the past 3 months
  df_min_p3m = lapply(seq_along(list_min_p3m), function(i) {
    list_min_p3m[[i]] %>% 
      dplyr::filter(date >= start_date) %>% 
      dplyr::select(one_of("min_p3m"))
  }) %>% 
    do.call("cbind",.) %>% 
    `colnames<-`(stations_dd_id)
  # minimum of the past 12 months
  df_min_p12m = lapply(seq_along(list_min_p12m), function(i) {
    list_min_p12m[[i]] %>% 
      dplyr::filter(date >= start_date) %>% 
      dplyr::select(one_of("min_p12m"))
  }) %>% 
    do.call("cbind",.) %>% 
    `colnames<-`(stations_dd_id)
  
  # mean of the past 3 months
  df_mean_p3m = lapply(seq_along(list_mean_p3m), function(i) {
    list_mean_p3m[[i]] %>% 
      dplyr::filter(date >= start_date) %>% 
      dplyr::select(one_of("mean_p3m"))
  }) %>% 
    do.call("cbind",.) %>% 
    `colnames<-`(stations_dd_id)
  # mean of the past 12 months
  df_mean_p12m = lapply(seq_along(list_mean_p12m), function(i) {
    list_mean_p12m[[i]] %>% 
      dplyr::filter(date >= start_date) %>% 
      dplyr::select(one_of("mean_p12m"))
  }) %>% 
    do.call("cbind",.) %>% 
    `colnames<-`(stations_dd_id)
  
  cat("the convertion of standard deviation and coffecient of varation
      have been left.")
  # sd
  df_sd = lapply(seq_along(sd_list), function(i) {
    sd_list[[i]] %>% 
      dplyr::filter(date >= start_date) %>% 
      dplyr::select(one_of("sd"))
  }) %>% 
    do.call("cbind",.) %>% 
    `colnames<-`(stations_dd_id)
  # cv
  df_cv = lapply(seq_along(cv_list), function(i) {
    cv_list[[i]] %>% 
      dplyr::filter(date >= start_date) %>% 
      dplyr::select(one_of("cv"))
  }) %>% 
    do.call("cbind",.) %>% 
    `colnames<-`(stations_dd_id)
  
  # combine all the predictors into a final list as output of the func.
  out <- list(min_p3m = df_min_p3m, mean_p3m = df_mean_p3m,
              min_p12m = df_min_p12m, mean_p12m = df_mean_p12m,
              sd = df_sd, cv = df_cv)
  
  return(out)
}
=======
  # find the stations with records for at least 36 months
  station_less36_mon <- select_min_fullmonths(out, min_months = 36)
  
  function_output <- list(
    output_ts = station_less36_mon,
    gauge_ids = station_less36_mon[, unique(gaugeid)]
  )
  
  return(function_output)
}

remove_records <- function(in_daily_paths) {
  #rbind everything
  
  #Remove all gauges with 0 values that have at least 99% of integer values as not reliable (see GRDC_NO 6140700 as example)
  GRDCtoremove_allinteger <- data.table(
    GRDC_NO = GRDCstatsdt[integerperc_o1800 >= 0.95 &
                            intermittent_o1800 == 1, GRDC_NO],
    flag = 'removed',
    comment = 'All integer discharge values'
  )
  
  #------ Remove stations based on examination of plots and data series
  daily_q_records_to_remove <- list(
    #C(gaugeid, 'removed' or 'to inspect' or'inspected', comment)
  ) %>%
    rbindlist %>%
    as.data.table %>%
    setnames(c('gaugeid', 'flag', 'comment'))
  
  
  gsim_q_records_to_remove <- list(
    #C(gaugeid, 'removed' or 'to inspect' or'inspected', comment)
  )%>%
    rbindlist %>%
    as.data.table %>%
    setnames(c('gaugeid', 'flag', 'comment'))
}

>>>>>>> origin/dev_merged

