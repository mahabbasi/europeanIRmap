#---------------------------- UTILITY FUNCTIONS --------------------------------
grdc_txt2csv <- function(path_list){
  
  start_date <<- start_date
  end_date <<- end_date
  
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
      .[dates >= start_date & dates <= end_date,] %>%
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
    st_set_crs(4236) %>% 
    st_transform(st_crs(4236)) %>%
    sf::st_intersection(., eu_countries)
  
  eu_final_sf <- joined_sf %>% 
    dplyr::select(-c(23:79)) 
  
}

select_dataset <- function(shp_path, dataset_name = NULL,
                           col_names = NULL){
  
  station_sf <- sf::st_read(dsn = shp_path, quiet = TRUE)
  
  interest_col = station_sf %>% 
    dplyr::filter(gaugdts == dataset_name) %>% 
    dplyr::pull(col_names)
  
  gauge_d =  station_sf %>% 
    dplyr::filter(gaugdts == dataset_name) %>% 
    dplyr::pull("gauge_d")
  
  return(list(interest_col, gauge_d))
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

select_spanish_stations <- function(in_paths, in_metadata) {
  
  start_date <<- start_date
  end_date <<- end_date
  
  spain_stations_id <- select_dataset(shp_path = shp_path,
                                             dataset_name = "spanish",
                                             col_names = "gauge_d")[[1]]
  
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
      .[dates >= start_date & dates <= end_date,]
    
    return(all_dt)
  }) %>% 
    rbindlist
  
  station_less36_mon <- select_min_fullmonths(out, min_months = 36) %>% 
    .[gaugeid %in% spain_stations_id,]
  # add the function to select spanish stations later
  
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
                                             col_names = "gauge_d")[[1]]
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
  
  start_date <<- start_date
  
  grdc_dd_id_old <- sf::st_read(shp_path, quiet = TRUE) %>% 
    dplyr::filter(gauge_d %in% old_grdc_stations_id) %>% 
    dplyr::pull(dd_id)
  
  grdc_new_names = c("dates", old_grdc_stations_id)
  
  # read the stations and produce a long format data table
  output <- data.table::fread(path) %>% 
    dplyr::select(one_of(c("date", grdc_dd_id_old))) %>%
    dplyr::mutate(date = as.Date(date, format = "%m/%d/%Y")) %>%
    dplyr::filter(date >= start_date) %>% 
    dplyr::rename_with(~grdc_new_names, .cols = everything()) %>%
    melt(id.vars='dates', variable.name='gaugeid')
  
  # find the stations with records for at least 36 months
  station_less36_mon <- select_min_fullmonths(output, min_months = 36)
  
  function_output <- list(
    output_ts = station_less36_mon,
    gauge_ids = station_less36_mon[, unique(gaugeid)]
  )
  
  return(function_output)
}

select_gsim_stations <- function(path, shp_path, dataset_name = "GSIM"){
  
  start_date <<- start_date
  
  # find the GSIM stations within the shapefile
  gsim_stations <- select_dataset(shp_path = shp_path,
                                  dataset_name = dataset_name,
                                  col_names = "dd_id")
  
  gsim_stations_dd_id = gsim_stations[[1]]
  gsim_stations_gauge_d = gsim_stations[[2]]
  
  gsim_new_names = c("dates", gsim_stations_gauge_d)
  
  output <- data.table::fread(path) %>% 
    dplyr::select(one_of(c("date", gsim_stations_dd_id))) %>%
    dplyr::mutate(date = as.Date(date, format = "%m/%d/%Y")) %>%
    dplyr::filter(date >= start_date) %>% 
    dplyr::rename_with(~gsim_new_names, .cols = everything()) %>%
    melt(id.vars='dates', variable.name='gaugeid')
  
  # find the stations with records for at least 36 months
  station_less36_mon <- select_min_fullmonths(output, min_months = 36)
  
  function_output <- list(
    output_ts = station_less36_mon,
    gauge_ids = station_less36_mon[, unique(gaugeid)]
  )
  
  return(function_output)
}

select_smires_stations <- function(path, shp_path){
  
  start_date <<- start_date
  end_date <<- end_date
  
  smires_stations_ls <- select_dataset(shp_path = shp_path,
                                       dataset_name = "smires",
                                       col_names = "dd_id")
  
  smires_stations_dd_id = smires_stations_ls[[1]]
  smires_stations_gauge_d = smires_stations_ls[[2]]
  
  out <- lapply(seq_along(path), function(i){
    
    #extract smires unique ID by formatting path
    gaugeno <- strsplit(basename(path[i]), '[.]')[[1]][1]
    index = which(smires_stations_dd_id == gaugeno )
    gaugeid <- smires_stations_gauge_d[index]
    
    
    # read the smires csv files for each station
    m <- data.table::fread(path[i], header = T, sep=",") %>%
      setnames('time', 'dates') %>% 
      setnames("streamflow", 'value') %>% 
      mutate(dates = as.Date(dates))
    
    # join the values of station by matching dates
    all_dt <- data.table(dates = seq.Date(
      as.Date(paste0(format(min(m$dates, na.rm=T), '%Y'),'-01-01')),
      as.Date(paste0(format(max(m$dates, na.rm=T), '%Y'),'-12-31')),
      "day"))
    
    all_dt <- merge(all_dt, m, by="dates", all.x=T) %>%
      .[dates >= start_date & dates <= end_date,] %>%
      .[, gaugeid := gaugeid] %>% 
      .[!is.na(gaugeid)]
    
    return(all_dt)
  }) %>% 
    rbindlist
  
  out[value %in% c(-999, -99, -9999, 999, 9999), value := NA] 
  
  # find the stations with records for at least 36 months
  station_less36_mon <- select_min_fullmonths(out, min_months = 36)
  
  function_output <- list(
    output_ts = station_less36_mon,
    gauge_ids = station_less36_mon[, unique(gaugeid)]
  )
  
  return(function_output)
}

select_corsica_stations <- function(path, shp_path){
  
  start_date <<- start_date
  end_date <<- end_date
  
  # find the corsica stations within the shapefile
  # corsica_stations_id <- select_dataset(shp_path = shp_path,
  #                                       dataset_name = "Corsica")[['gauge_d']]
  
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
    
    all_dt <- merge(all_dt, m, by="dates", all.x=T) %>%
      .[dates >= start_date & dates <= end_date,] %>%
      .[, gaugeid := gaugeid]
    
    return(all_dt)
  })  %>% 
    rbindlist
  
  out[value %in% c(-999, -99, -9999, 999, 9999), value := NA] 
  out[, value := value/1000]

  # find the stations with records for at least 36 months
  station_less36_mon <- select_min_fullmonths(out, min_months = 36)
  
  function_output <- list(
    output_ts = station_less36_mon,
    gauge_ids = station_less36_mon[, unique(gaugeid)]
  )
  
  return(function_output)
}

select_italian_emr_stations <- function(path, shp_path){
  
  start_date <<- start_date
  end_date <<- end_date
  
  # find the emr stations within the shapefile
  emr_stations_ls <- select_dataset(shp_path = shp_path,
                                    dataset_name = "arpae",
                                    col_names = "name")
  
  emr_stations_id = emr_stations_ls[[1]]
  emr_stations_gauge_d = emr_stations_ls[[2]]
  
  emr_stations_id <- emr_stations_id[! emr_stations_id %in% "emr:Salsominore"]
  emr_stations_id <- c(emr_stations_id, "emr:Salsominiore")
  
  out <- lapply(seq_along(path), function(i){
    
    gaugeno <- gsub(".csv", "", basename(path[i])) %>%
      gsub("emr_", "emr:", .) %>% 
      gsub("_", " ", .)
    
    index = which(emr_stations_id == gaugeno)
    gaugeid <- emr_stations_gauge_d[index]
    
    m <- fread(path[i], header = T, skip = 14) %>% .[,-2] %>% 
      rename(dates = "Inizio validità (UTC)")  %>% 
      setnames("Portata media giornaliera (M**3/S)", 'value') %>% 
      dplyr::mutate(dates = as.Date(dates))
    
    # join the values of station by matching dates
    all_dt <- data.table(dates = seq.Date(
      as.Date(paste0(format(min(m$dates, na.rm=T), '%Y'),'-01-01')),
      as.Date(paste0(format(max(m$dates, na.rm=T), '%Y'),'-12-31')),
      "day"))
    
    all_dt <- merge(all_dt, m, by="dates", all.x=T) %>%
      .[dates >= start_date & dates <= end_date,] %>%
      .[, gaugeid := gaugeid] %>% 
      .[!is.na(gaugeid)]
    
    return(all_dt)
  }) %>% 
    rbindlist
  
  # find the stations with records for at least 36 months
  station_less36_mon <- select_min_fullmonths(out, min_months = 36)
  
  function_output <- list(
    output_ts = station_less36_mon,
    gauge_ids = station_less36_mon[, unique(gaugeid)]
  )
  
  return(function_output)
}

select_italian_ispra_stations <- function(path, shp_path){
  
  start_date <<- start_date
  end_date <<- end_date
  
  # find the corsica stations within the shapefile
  ispra_stations_ls <- select_dataset(shp_path = shp_path,
                                      dataset_name = "ispra",
                                      col_names = "name")
  
  ispra_stations_id = ispra_stations_ls[[1]]
  ispra_stations_gauge_d = ispra_stations_ls[[2]]
  
  out <- lapply(seq_along(path), function(i){
    # find the station name from the path
    gaugeno <- gsub(".csv", "", basename(path[i])) %>%
      gsub("_", ":", .)
    
    index = which(ispra_stations_id == gaugeno)
    gaugeid <- ispra_stations_gauge_d[index]
    
    m <- fread(path[i], header = T, skip = 15) %>% 
      rename(dates = "flowvalue") %>% 
      setnames("dateTime", 'value')
    
    # join the values of station by matching dates
    all_dt <- data.table(dates = seq.Date(
      as.Date(paste0(format(min(m$dates, na.rm=T), '%Y'),'-01-01')),
      as.Date(paste0(format(max(m$dates, na.rm=T), '%Y'),'-12-31')),
      "day"))
    
    all_dt <- merge(all_dt, m, by="dates", all.x=T) %>%
      .[dates >= start_date & dates <= end_date,] %>%
      .[, gaugeid := gaugeid]
    
    return(all_dt)
  }) %>% 
    rbindlist
  
  # find the stations with records for at least 36 months
  station_less36_mon <- select_min_fullmonths(out, min_months = 36) %>% 
    .[gaugeid %in% ispra_stations_gauge_d,]
  
  function_output <- list(
    output_ts = station_less36_mon,
    gauge_ids = station_less36_mon[, unique(gaugeid)]
  )
  
  return(function_output)
}

select_italian_arpal_stations <- function(path, shp_path){
  
  start_date <<- start_date
  end_date <<- end_date
  
  # find the arpal stations within the shapefile
  arpal_stations_ls <- select_dataset(shp_path = shp_path,
                                      dataset_name = "arpal",
                                      col_names = "name")
  arpal_stations_id = arpal_stations_ls[[1]]
  arpal_stations_gauge_d = arpal_stations_ls[[2]]
  
  out <- lapply(seq_along(path), function(i){
    #extract emr unique ID by formatting path
    gaugeno <- gsub(".csv", "", basename(path[i]))
    
    index = which(arpal_stations_id == gaugeno)
    gaugeid <- arpal_stations_gauge_d[index]

    m <- fread(path[i], header = TRUE) %>% .[, -1] %>% 
      mutate(date = as.Date(date, format = "%Y_%m_%d"),
             Q = suppressWarnings(as.numeric(Q))) %>% 
      setnames(c("date", "Q"),
               c("dates", 'value'))
    
    # join the values of station by matching dates
    all_dt <- data.table(dates = seq.Date(
      as.Date(paste0(format(min(m$dates, na.rm=T), '%Y'),'-01-01')),
      as.Date(paste0(format(max(m$dates, na.rm=T), '%Y'),'-12-31')),
      "day"))
    
    # join the values of station by matching dates
    all_dt <- merge(all_dt, m, by="dates", all.x=T) %>%
      .[dates >= start_date & dates <= end_date,] %>%
      .[, gaugeid := gaugeid] %>% 
      .[!is.na(gaugeid)]
    
    return(all_dt)
  }) %>% 
    rbindlist
  
  # find the stations with records for at least 36 months
  station_less36_mon <- select_min_fullmonths(out, min_months = 36)
  
  function_output <- list(
    output_ts = station_less36_mon,
    gauge_ids = station_less36_mon[, unique(gaugeid)]
  )
  
  return(function_output)
}

select_italian_arpas_stations <- function(path, shp_path){
  
  # call the start and end date of the period 
  start_date <<- start_date
  end_date <<- end_date
  # all_dt <- data.table(dates = seq.Date(as.Date("1901-01-01"),
  #                                       as.Date("2022-12-31"),
  #                                       "day"))
  
  # find the arpas stations within the shapefile
  arpas_stations_ls <- select_dataset(shp_path = shp_path,
                                      dataset_name = "arpas",
                                      col_names = "name")
  
  arpas_stations_id = arpas_stations_ls[[1]]
  arpas_stations_gauge_d = arpas_stations_ls[[2]]
  
  out <- lapply(seq_along(path), function(i){
    
    # find the station name from the path
    gaugeno <- gsub(".csv", "", basename(path[i])) %>%
      gsub("_",".", .) %>% 
      paste0("sar:", .)
    
    index = which(arpas_stations_id == gaugeno)
    gaugeid <- arpas_stations_gauge_d[index]
    
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
      .[dates >= start_date & dates <= end_date,] %>%
      .[, gaugeid := gaugeid] %>% 
      .[!is.na(gaugeid)]
    
    return(all_dt)
  }) %>% 
    rbindlist
  
  # find the stations with records for at least 36 months
  station_less36_mon <- select_min_fullmonths(out, min_months = 36)
  
  function_output <- list(
    output_ts = station_less36_mon,
    gauge_ids = station_less36_mon[, unique(gaugeid)]
  )
  
  return(function_output)
}

remove_records <- function(in_dt_daily) {
  #rbind everything and make a copy of inputs dts
  if (inherits(in_dt_daily, "list")) {
    gaugetab = in_dt_daily %>% 
      rbindlist() %>% 
      .[,gaugeid := as.character(gaugeid)] %>% 
      copy(.)
  } else {
    gaugetab = in_dt_daily %>% 
      .[,gaugeid := as.character(gaugeid)] %>% 
      copy(.)
  }

  
  #------ Remove stations based on examination of plots and data series
  
  daily_q_records_to_remove = list(
    # c("station_id", "dataset_name","comments", "flag", "start_date", "end_date")
    c("UA_0000023", "GSIM","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("sm_1281", "smires","abrupt decrease to zero", 'months with zeros', "2009-01-01", "2013-12-31"),
    c("sm_1245", "smires","abrupt decrease to zero", 'months with zeros','1992-01-01', '1992-12-31'),
    c("sm_1181", "smires", "getting start with zero and then never been again", 'remove entire record', NA, NA),
    c("sm_1138", "smires", "getting start with zero and then never been again", 'remove entire record', NA, NA),
    c("sm_1124", "smires", "abrupt decrease to zero", 'months with zeros', NA, NA),
    c("sm_1118", "smires", "abrupt decrease to zero", 'months with zeros', '1987-04-01', '1987-04-30'),
    c("sm_1080", "smires", "abrupt decrease to zero", 'remove entire record', NA, NA),
    c("sm_1048", "smires", "abrupt decrease to zero", 'months with zeros', '1993-10-01', '1993-10-31'),
    c("sm_1038", "smires", "abrupt decrease to zero", 'months with zeros', '1984-09-01', '1984-09-30'),
    c("arpal_2004", "arpal","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("arpal_2002", "arpal","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("6854103", "GRDC","abrupt decrease to zero & round records", 'remove entire record', NA, NA),
    c("6642100", "GRDC","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("6503352", "GRDC","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("6444500", "GRDC","abrupt decrease to zero & round records", 'remove entire record', NA, NA),
    c("6337530", "GRDC","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("6233410", "GRDC","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("6125310", 'GRDC', 'abrupt decrease to anomalously low values and zeros at the end of the series', 'months with zeros', '2009-01-01', NA),
    c("6123760", 'GRDC', 'abrupt decrease to anomalously low values and zeros at the end of the series', 'months with zeros', '2009-01-01', NA),
    c("6123641", "GRDC","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("6123370", 'GRDC', 'abrupt decrease to anomalously low values and zeros at the end of the series', 'months with zeros', '2009-01-01', NA),
    c("9263", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("9231", "spanish_stations","start with zero values but never exprience again", 'remove entire record', NA, NA),
    c("9181", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("9170", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("9157", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("9135", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("9096", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("9090", "spanish_stations","abrupt decrease to zero and rounded records to two digits", 'remove entire record', NA, NA),
    c("9005", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("8132", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("8130", "spanish_stations","start with zero values but never exprience again", 'remove entire record', NA, NA),
    c("8071", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("7102", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("7062", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("7050", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("7006", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("7004", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("7003", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("5084", "spanish_stations","abrupt decrease to zero", 'months with zeros', '2004-01-01', '2004-03-31'),
    c("5043", "spanish_stations","abrupt decrease to zero and rounded records to two digits", 'remove entire record', NA, NA),
    c("5039", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("5003", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("5002", "spanish_stations","abrupt decrease to zero and rounded records to two digits", 'remove entire record', NA, NA),
    c("4122", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("3285", "spanish_stations","abrupt decrease to zero and rounded records to two digits", 'remove entire record', NA, NA),
    c("3283", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("3270", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("3263", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("3252", "spanish_stations","abrupt decrease to zero and duplicated values in months", 'remove entire record', NA, NA),
    c("3155", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("3154", "spanish_stations","records is too incomplete and unreliable", 'remove entire record', NA, NA),
    c("3113", "spanish_stations","records is too incomplete and unreliable", 'remove entire record', NA, NA),
    c("2125", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("2064", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("1544", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("1464", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("1196", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("6", "rbis","abrupt decrease to zero", 'months with zeros', NA, NA),
    c('6854510', 'GRDC', "abrupt decrease to zero. rounded to nearest m3", 'remove entire record', NA, NA),
    c('6683300', 'GRDC', "rounded to nearest m3", 'remove entire record', NA, NA),
    c('6683200', 'GRDC', "abrupt decrease to zero. rounded to nearest m3", 'remove entire record', NA, NA),
    c('6683010', 'GRDC', "abrupt decrease to zero. rounded to nearest m3", 'remove entire record', NA, NA),
    c('6682300', 'GRDC', "abrupt decrease to zero. rounded to nearest m3", 'remove entire record', NA, NA),
    c('6444400', 'GRDC', "abrupt decrease to zero. rounded to nearest m3", 'remove entire record', NA, NA),
    c('6444250', 'GRDC', "abrupt decrease to zero. rounded to nearest m3", 'remove entire record', NA, NA),
    c('6373911', 'GRDC', "abrupt decrease to zero. records is too incomplete and  unreliable", 'remove entire record', NA, NA),
    c('6373460', 'GRDC', "abrupt decrease to zero. records is too incomplete and  unreliable", 'remove entire record', NA, NA),
    c('6373444', 'GRDC', "abrupt decrease to zero. records is too incomplete and  unreliable", 'remove entire record', NA, NA),
    c('6373224', 'GRDC', "abrupt decrease to zero. records is too incomplete and  unreliable", 'remove entire record', NA, NA),
    c('6233740', 'GRDC', "abrupt decrease to zero", 'remove entire record', NA, NA),
    c('6221660', 'GRDC', "an anomalously jump in the period", 'remove entire record', NA, NA),
    c('6125300', 'GRDC', "abrupt decrease to zero. rounded to nearest m3", 'remove entire record', NA, NA),
    c('9292', 'spanish_stations', "abrupt decrease to zero. rounded to nearest m3", 'remove entire record', NA, NA),
    c('9145', 'spanish_stations', "start with zero values but later an anomalously jump in the period", 'remove entire record', NA, NA),
    c('9143', 'spanish_stations', "an anomalously jump in the period. record is unreliable", 'remove entire record', NA, NA),
    c('9105', 'spanish_stations', "abrupt decrease to zero. rounded to two digits", 'remove entire record', NA, NA),
    c('7007', 'spanish_stations', "abrupt decrease to zero. rounded to two digits", 'remove entire record', NA, NA),
    c('3203', 'spanish_stations', "abrupt decrease to zero. rounded to two digits", 'remove entire record', NA, NA),
    c('3128', 'spanish_stations', "abrupt decrease to zero. rounded to two digits", 'remove entire record', NA, NA),
    c('3019', 'spanish_stations', "abrupt decrease to zero. rounded to two digits", 'remove entire record', NA, NA),
    c('3016', 'spanish_stations', "abrupt decrease to zero. records is too incomplete and  unreliable", 'remove entire record', NA, NA),
    c(5, 'rbis', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(1264, 'spanish_stations','abrupt decrease to zero', 'months with zeros', '2007-01-01', '2007-01-31'),
    c(1295, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(1398, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero', 'months with zeros', NA, NA),
    c(1431, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(1438, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(1446, 'spanish_stations', 'zero flows may be due to values being rounded to two decimals',  'months with zeros', NA, NA),
    c(1485, 'spanish_stations', 'zero flows may be due to values being rounded to two decimals. odd patterns probably due to regulation',  'months with zeros', NA, NA),
    c(1519, 'spanish_stations', 'zero flows may be due to values being rounded to two decimals and rating curve discontinuity',  'months with zeros', NA, NA),
    c(1520, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(1542, 'spanish_stations', 'zero flows may be due to values being rounded to two decimals and rating curve discontinuity',  'months with zeros', NA, NA),
    c(2000, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero','months with zeros', NA, '2001-01-01'),
    c(2003, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero', 'remove entire record', NA, NA),
    c(2006, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero', 'months with zeros', NA, NA),
    c(2009, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero', 'months with zeros', NA, '2001-01-01'),
    c(2015, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero', 'months with zeros', NA, NA),
    c(2016, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero', 'months with zeros', NA, '1983-01-01'),
    c(2018, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(2030, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(2031, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, '1985-01-01'),
    c(2033, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero', 'remove entire record', NA, NA),
    c(2035, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero', 'months with zeros', '2000-08-01', '2000-09-30'),
    c(2040, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero', 'months with zeros', NA, NA),
    c(2053, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero', 'months with zeros', NA, '2001-01-01'),
    c(2054, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, '1986-01-01'),
    c(2056, 'spanish_stations', 'zero flows may be due to values being rounded to two decimals at the beginning of the time series','months with zeros', NA, '1986-01-01'),
    c(2066, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(2068, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(2081, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero','months with zeros', NA, '2008-01-01'),
    c(2083, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero', 'months with zeros', NA, NA),
    c(2093, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(2097, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero', 'months with zeros', NA, NA),
    c(2102, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero', 'months with zeros', NA, NA),
    c(2109, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(2129, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero', 'months with zeros', NA, '2001-01-01'),
    c(2144, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(2145, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(2149, 'spanish_stations', 'zero flows may be due to values being rounded to two decimals',  'months with zeros', NA, NA),
    c(3001, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(3009, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(3015, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(3030, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(3065, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero', 'months with zeros', NA, NA),
    c(3067, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero', 'months with zeros', NA, NA),
    c(3127, 'spanish_stations', 'zero flows may be due to values being rounded to one decimals',  'remove entire record', NA, NA),
    c(3142, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(3144, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros','2008-12-01', '2009-12-31'),
    c(3146, 'spanish_stations', 'zero flows may be due to values being rounded to two decimals', 'months with zeros', NA, '2010-01-01'),
    c(3152, 'spanish_stations', 'rounded to nearest m3 or first decimal', 'remove entire record', NA, NA),
    c(3159, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero', 'months with zeros', NA, NA),
    c(3162, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero', 'months with zeros', NA, NA),
    c(3175, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(3180, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero','months with zeros', NA, '2010-01-01'),
    c(3189, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(3190, 'spanish_stations', 'zero flows may be due to values being rounded to two decimals', 'remove entire record', NA, NA),
    c(3195, 'spanish_stations', 'zero flows may be due to values being rounded to two decimals', 'remove entire record', NA, NA),
    c(3196, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(3200, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', '2005-01-01', '2005-10-31'),
    c(3212, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, '1987-01-01'),
    c(3262, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero', 'remove entire record', NA, NA),
    c(3269, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero', 'remove entire record', NA, NA),
    c(3940, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(4009, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero', 'remove entire record', NA, NA),
    c(4013, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero', 'remove entire record', NA, NA),
    c(4014, 'spanish_stations', 'record is too incomplete', 'remove entire record', NA, NA),
    c(4105, 'spanish_stations', 'record is too incomplete', 'remove entire record', NA, NA),
    c(4131, 'spanish_stations', 'record is too incomplete', 'remove entire record', NA, NA),
    c(4140, 'spanish_stations', 'record is too incomplete', 'remove entire record', NA, NA),
    c(4164, 'spanish_stations', 'record is too incomplete', 'remove entire record', NA, NA),
    c(4165, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, '2004-01-01'),
    c(4229, 'spanish_stations', 'record is too incomplete', 'remove entire record', NA, NA),
    c(4248, 'spanish_stations', 'record is too incomplete', 'remove entire record', NA, NA),
    c(4283, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(4287, 'spanish_stations', 'record is too incomplete', 'remove entire record', NA, NA),
    c(5008, 'spanish_stations', 'zero flows may be due to values being rounded to two decimals', 'remove entire record', NA, NA),
    c(5016, 'spanish_stations', 'record is too incomplete', 'months with zeros', NA, '2015-01-01'),
    c(5023, 'spanish_stations', 'record is too incomplete', 'months with zeros', NA, '2015-01-01'),
    c(5028, 'spanish_stations', 'record is too incomplete', 'remove entire record', NA, NA),
    c(5047, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero', 'months with zeros', NA, NA),
    c(5082, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero', 'months with zeros', NA, '1985-01-01'),
    c(5083, 'spanish_stations', 'abrupt decrease and increase', 'months with zeros', NA, NA),
    c(5086, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(5101, 'spanish_stations', 'record is too incomplete', 'months with zeros', '2010-01-01', '2010-12-31'),
    c(5133, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero', 'months with zeros', NA, '1983-01-01'),
    c(5138, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(7043, 'spanish_stations', 'rating curve discontinuity', 'remove entire record', NA, NA),
    c(7121, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(8028, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(8032, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(8042, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(8074, 'spanish_stations', 'rating curve discontinuity', 'months with zeros', NA, '1988-01-01'),
    c(8148, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(9014, 'spanish_stations', 'rating curve discontinuity','months with zeros', NA, NA),
    c(9018, 'spanish_stations', 'rating curve discontinuity','months with zeros', NA, NA),
    c(9035, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(9039, 'spanish_stations', 'rating curve discontinuity', 'months with zeros', NA, '2010-01-01'),
    c(9049, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', '1986-09-01','1986-09-30'),
    c(9057, 'spanish_stations', 'rating curve discontinuity', 'months with zeros', NA, '2011-01-01'),
    c(9088, 'spanish_stations', 'rating curve discontinuity', 'months with zeros', NA, NA),
    c(9095, 'spanish_stations', 'rating curve discontinuity', 'months with zeros', NA, NA),
    c(9097, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(9109, 'spanish_stations', 'rating curve discontinuity', 'remove entire record', NA, NA),
    c(9115, 'spanish_stations', 'rating curve discontinuity', 'months with zeros', NA, NA),
    c(9118, 'spanish_stations', 'rating curve discontinuity', 'months with zeros', NA, '2010-01-01'),
    c(9124, 'spanish_stations', 'rating curve discontinuity', 'months with zeros', NA, '2010-01-01'),
    c(9138, 'spanish_stations', 'rating curve discontinuity', 'remove entire record', NA, NA),
    c(9144, 'spanish_stations', 'rating curve discontinuity', 'remove entire record', NA, NA),
    c(9159, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', '1987-01-01', '1996-01-01'),
    c(9182, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(9186, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(9201, 'spanish_stations', 'rating curve discontinuity','months with zeros', NA, NA),
    c(9257, 'spanish_stations', 'rating curve discontinuity','remove entire record', NA, NA),
    c(9259, 'spanish_stations', 'rating curve discontinuity','months with zeros', NA, NA),
    c(9287, 'spanish_stations', 'rating curve discontinuity','months with zeros', NA, NA),
    c(9321, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(6114500, 'GRDC',  'zero flows may be due to values being rounded to two decimals',  'months with zeros', NA, NA),
    c(6119100, 'GRDC', 'zero flows may be due to values being rounded to two decimals', 'months with zeros', NA, NA),
    c(6123461, 'GRDC', 'abrupt decrease to anomalously low values and zeros at the end of the series', 'months with zeros', '2009-01-01', NA),
    c(6123501, 'GRDC', 'abrupt decrease to anomalously low values and zeros at the end of the series', 'months with zeros', '2009-01-01', NA),
    c(6124300, 'GRDC', 'abrupt decrease to anomalously low values and zeros at the end of the series', 'months with zeros', '2009-01-01', NA),
    c(6124430, 'GRDC', 'abrupt decrease to anomalously low values and zeros at the end of the series', 'months with zeros', '2009-01-01', NA),
    c(6124440, 'GRDC', 'abrupt decrease to anomalously low values and zeros at the end of the series', 'months with zeros', '2009-01-01', NA),
    c(6124502, 'GRDC', 'abrupt decrease to anomalously low values and zeros at the end of the series', 'months with zeros', '2009-01-01', NA),
    c(6125360, 'GRDC', 'abrupt decreases to zero at the end of the series', 'months with zeros', '1996-01-01', NA),
    c(6125440, 'GRDC', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(6128050, 'GRDC', 'zero flows may be due to values being rounded to two decimals. flow regulation started within record',  'months with zeros', '1989-01-01', NA),
    c(6128101, 'GRDC', 'abrupt decrease to anomalously low values and zeros at the end of the series', 'months with zeros', '2009-01-01', NA),
    c(6128220, 'GRDC', 'abrupt decrease to zeor after gap in series', 'months with zeros', '2007-04-01', '2007-05-01'),
    c(6139340, 'GRDC', 'rounded to nearest m3', 'remove entire record', NA, NA),
    c(6139501, 'GRDC', 'abrupt decrease to anomalously low values and zeros at the end of the series', 'months with zeros', '2009-01-01', NA),
    c(6139700, 'GRDC', 'rounded to nearest m3', 'remove entire record', NA, NA),
    c(6139850, 'GRDC', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(6144490, 'GRDC', 'rounded to nearest m3', 'remove entire record', NA, NA),
    c(6172028, 'GRDC', 'decrease to low discharge may be plausible but zero flows may be due to values being rounded to two decimals', 'months with zeros', NA, NA),
    c(6172031, 'GRDC', 'zero flows may be due to values being rounded to two decimals',  'remove entire record', NA, NA),
    c(6233128, 'GRDC', 'low-flows could be miscalibrated zeros', 'remove entire record', NA, NA),
    c(6233250, 'GRDC', 'rounded to nearest m3', 'remove entire record', NA, NA),
    c(6233414, 'GRDC', 'rounded to nearest m3', 'remove entire record', NA, NA),
    c(6242810, 'GRDC', 'abrupt decrease to zero. probably due to regulation for channel maintenance or other', 'months with zeros', NA, NA),
    c(6273752, 'GRDC', 'zero flows may be due to values being rounded to two decimal', 'months with zeros', NA, NA),
    c(6273900, 'GRDC', 'rounded to nearest m3', 'remove entire record', NA, NA),
    c(6274550, 'GRDC', 'rounded to nearest m3', 'remove entire record', NA, NA),
    c(6274655, 'GRDC', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(6373040, 'GRDC', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(6373050, 'GRDC', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(6373200, 'GRDC', 'abrupt decrease to zero. unreliable record', 'remove entire record', NA, NA),
    c(6373215, 'GRDC', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(6373217, 'GRDC', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(6373306, 'GRDC', 'rounded to nearest m3', 'remove entire record', NA, NA),
    c(6442300, 'GRDC', 'rounded to nearest m3 until 1992', 'months with zeros', NA, '1992-01-01'),
    c(6444350, 'GRDC', 'rounded to nearest m3 until 1992', 'months with zeros', NA, '1992-01-01'),
    c(6458420, 'GRDC', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(6574154, 'GRDC', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(6680400, 'GRDC', 'rounded to nearest m3', 'remove entire record', NA, NA),
    c(6680410, 'GRDC', 'rounded to nearest m3', 'remove entire record', NA, NA),
    c(6682500, 'GRDC', 'rounded to nearest m3', 'remove entire record', NA, NA),
    c(6729180, 'GRDC', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(6729225, 'GRDC', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(6729310, 'GRDC', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(6731556, 'GRDC', 'abrupt decrease to zero at the end of the series', 'months with zeros', '2015-01-01', NA),
    c(6731815, 'GRDC', 'abrupt decrease to zero at the end of the series', 'months with zeros', '2006-01-01', NA),
    c(6731940, 'GRDC', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(6854580, 'GRDC', 'rounded to one decimal until 1992', 'months with zeros', NA, '1990-01-01'),
    c(6854591, 'GRDC', 'rounded to nearest m3', 'remove entire record', NA, NA),
    c(6854660, 'GRDC', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(6854708, 'GRDC', 'rounded to nearest m3', 'remove entire record', NA, NA),
    c(6854711, 'GRDC', 'rounded to nearest m3', 'remove entire record', NA, NA),
    c(6854722, 'GRDC', 'rounded to first decimal', 'remove entire record', NA, NA),
    c(6871100, 'GRDC', 'rounded to first decimal', 'remove entire record', NA, NA),
    c(6731940, 'GRDC', 'abrupt decrease to zero', 'remove entire record', NA, NA),
    c(6401120, 'GRDC', 'anomalously increased during the period', 'remove entire record', NA, NA),
    c('arpae_1013', 'emr', 'zero flows may be due to values being rounded to two decimal or to regulation', 'keep', NA, NA),
    c('arpae_1016', 'emr', 'zero flows may be due to values being rounded to two decimal or to regulation', 'keep', NA, NA),
    c('arpae_1017', 'emr', 'zero flows may be due to values being rounded to two decimal or to regulation', 'keep', NA, NA),
    c('arpae_1018', 'emr', 'zero flows may be due to values being rounded to two decimal or to regulation', 'keep', NA, NA),
    c('arpae_1020', 'emr', 'zero flows may be due to values being rounded to two decimal or to regulation', 'keep', NA, NA),
    c('arpae_1024', 'emr', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c('arpae_1032', 'emr', 'zero flows may be due to values being rounded to two decimal or to regulation', 'keep', NA, NA),
    c('arpae_1034', 'emr', 'zero flows may be due to values being rounded to two decimal or to regulation', 'keep', NA, NA),
    c('arpae_1035', 'emr', 'zero flows may be due to values being rounded to two decimal or to regulation', 'keep', NA, NA),
    c('arpae_1057', 'emr', 'zero flows may be due to values being rounded to two decimal or to regulation', 'keep', NA, NA),
    c('arpae_1072', 'emr', 'zero flows may be due to values being rounded to two decimal or to regulation', 'keep', NA, NA),
    c('arpae_1073', 'emr', 'zero flows may be due to values being rounded to two decimal or to regulation', 'keep', NA, NA),
    c('arpae_1079', 'emr', 'zero flows may be due to values being rounded to two decimal or to regulation', 'keep', NA, NA),
    c('arpal_2010', 'arpal', 'zero flows may be due to values being rounded to two decimal or to regulation', 'months with zeros', NA, NA),
    c('arpal_2028', 'arpal', 'zero flows may be due to values being rounded to two decimal or to regulation', 'months with zeros', NA, NA),
    c('arpas_3007', 'arpas', 'zero flows may be due to values being rounded to two decimal or to regulation', 'months with zeros', NA, NA),
    c('HS103', 'rbis', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c('PT_0000025', 'GSIM', 'zero flows may be due to values being rounded to two decimal or to regulation', 'months with zeros', NA, NA),
    c('sm_1003', 'smires', 'low-flows should probably be zero flows starting in 2016', 'change to zero', '2016-08-24', NA),
    c('sm_1171', 'smires', 'record is too incomplete', 'months with zeros', '2007-05-31', NA),
    c('sm_1174', 'smires', 'record is too incomplete', 'months with zeros', '2011-01-01', NA),
    c('sm_1180', 'smires', 'record is too incomplete', 'months with zeros', '2011-01-01', NA)
  )  %>% 
    do.call('rbind', .) %>% 
    as.data.table() %>% 
    `colnames<-`(c("station_id", "dataset_name","comments", "flag", "start_date", "end_date"))
  
  # Add the year and month columns to the data table to expand NA for the whole month
  gaugetab[, year_month := format(dates, "%Y-%m")]
  # Determine the different situations of the flags to exculd the gauges  
  gaugeid_removed <- daily_q_records_to_remove[flag == "remove entire record", station_id]
  
  gaugeid_zero_all <- daily_q_records_to_remove[flag == "months with zeros" &
                                      is.na(start_date) & is.na(end_date), station_id]
  
  gaugeid_zero_part <- daily_q_records_to_remove[flag == "months with zeros" &
                                       !is.na(start_date) | !is.na(end_date)] %>% 
    .[is.na(end_date), end_date := '2019-12-31'] %>% 
    .[is.na(start_date), start_date := '1981-01-01'] %>% 
    .[, c('start_date', 'end_date') := lapply(.SD, as.Date),
      .SDcols = c('start_date', 'end_date')]
  
  gaugeid_lowflow_zero <- daily_q_records_to_remove[flag == "change to zero" &
                                                      (!is.na(start_date) | !is.na(end_date))] %>% 
    .[is.na(end_date), end_date := '2019-12-31'] %>% 
    .[is.na(start_date), start_date := '1981-01-01'] %>% 
    .[, c('start_date', 'end_date') := lapply(.SD, as.Date),
      .SDcols = c('start_date', 'end_date')]
  
  # Implement the flags to the data table
  gaugetab[gaugeid %in% gaugeid_removed, value :=NA]
  
  gaugetab[gaugeid %in% gaugeid_zero_all & value == 0, value :=NA]
  
  for (iterator in 1:nrow(gaugeid_lowflow_zero)) {
    
    gaugeid_item <- gaugeid_lowflow_zero[iterator, station_id]
    start_date_item <- gaugeid_lowflow_zero[iterator, start_date]
    end_date_item <- gaugeid_lowflow_zero[iterator, end_date]
    gaugetab[gaugeid %in% gaugeid_item & 
               dates >= start_date_item & 
               dates <= end_date_item &
               value < 0.004,
             value:=0]
    
  }

  # Iterate over the data table with the gauges that modification requried in the middle
  # the period
  
  for (iterator in 1:nrow(gaugeid_zero_part)) {
    
    gaugeid_item <- gaugeid_zero_part[iterator, station_id]
    start_date_item <- gaugeid_zero_part[iterator, start_date]
    end_date_item <- gaugeid_zero_part[iterator, end_date]
    gaugetab[gaugeid %in% gaugeid_item & 
               dates >= start_date_item & 
               dates <= end_date_item &
               value == 0,
             value:=NA]
    
  }
  
  # Define the days of NA grouped by gauges and year_month then,
  # make all of the months NA with at least on one NA day
  gaugetab[, value := if (any(is.na(value))) NA else value, by = .(gaugeid, year_month)]
  

  
  # Select the non-NA records and return it
  output <- gaugetab[!is.na(value)] %>%
    .[,c('gaugeid', 'dates', 'value')]
  
  function_output <- list(output_ts = output,
                          gaugeids_removed = gaugeid_removed)
  
  return(function_output)

}

split_file_random = function(source_path, set1_path, set2_path, overlap_percentage = 0.1){
  
  # Get the files in the source path
  files_list = list.files(path = source_path, 
                          pattern = "*.png",
                          full.names = TRUE)
  
  # Set a seed for reproducibility
  set.seed(2023)
  overlap_size = round(length(files_list) * overlap_percentage)
  
  # Assign files randomly to two sets
  ovelap_indices = sample(length(files_list), size = overlap_size)
  sets_indices = sample(length(files_list), size = length(files_list)/2)
  
  # Get the files and unique them
  set1_files = files_list[unique(c(sets_indices, ovelap_indices))]
  set2_files =  files_list[-sets_indices]
  set2_files_ovelap = files_list[ovelap_indices]
  
  # Unique the paths
  set2_files_unique = unique(c(set2_files, set2_files_ovelap))
  
  # Copy files to the subfolders
  file.copy(set1_files, set1_path)
  file.copy(set2_files_unique, set2_path)
  
}

final_gauge_shp <- function(shp_path, gaugeids_removed){
  
  shapefiles <- sf::st_read(dsn = shp_path, quiet = TRUE)
  
  shp_filtered <- shapefiles %>% 
    dplyr::filter(!gauge_d %in% gaugeids_removed)
  
  return(shp_filtered)
}

# -------------- Extract predictors, HR, LR, Statics ------------------------------------
### High Resolution predictors -----
extract_HR_streamflow = function(path, gauges_shp, start_year = 1980,
                                 end_year = 2019, out_dir = NULL){
  
  years = seq(start_year, end_year)
  
  filename = vector("character")
  for (i in 1:length(years)) {
    
    for (j in 1:12) {
      filename[(i-1)*12 + j] = paste0("15sec_dis_", j, "_", years[i], ".tiff")
    }
  }
  
  
  WaterGap_mat = matrix(NA, nrow = length(filename),
                        ncol = dim(gauges_shp)[1])
  for (i in seq_along(filename)) {
    datapath = file.path(path, filename[i])
    r = raster::raster(datapath)
    rasValue = raster::extract(r, gauges_shp)
    WaterGap_mat[i,] = rasValue
    cat("The extraction of the streamflow for the gauging stations in",
        filename[i], "file was done.\n")
  }
  
  colnames(WaterGap_mat) = gauges_shp$gauge_d
  
  WaterGap_mat %>%
    data.table::as.data.table() %>%
    data.table::fwrite(., file = file.path(out_dir,
                                             "waterGAP_streamflow_stations.csv"))
  
  return(WaterGap_mat)
}

#' add_date: this function is used to generate a sequence monthly date within our time span.
#' 
#' @param mat (matrix) the time series of WaterGAP streamflow for all the stations.
#'  
add_date <- function(mat){
  mat %>% 
    as.data.frame(.) %>% 
    replace(is.na(.), 0) %>%
    mutate(date = seq.Date(as.Date("1980-01-01"),
                           as.Date("2019-12-31"),
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
sd_mon <- function(tbl, start_date = "1980-01-01"){
  tbl %>%
    filter(date >= as.Date(start_date)) %>%
    as.vector(.) %>% .$Q %>% matrix(., ncol = 12, byrow = TRUE) %>%
    as.data.frame(.) %>% `colnames<-`(month.abb) %>%
    summarise_all(., sd)
} 

#' calculating the average of streamflow for 12 months for the gauging stations.
#' 
#'  @param tbl (vector) the gauging station streamflow vector in either dataframe or tibble class. 
#'  @param start_date (character) the start date of the time period
#'
mean_mon <- function(tbl, start_date = "1980-01-01"){
  tbl %>%
    filter(date >= as.Date(start_date)) %>%
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
compute_highres_predictors <- function(watergap_raw_path){
  
  start_date <<- start_date
  end_date <<- end_date
  # read the monthly downscaled watergap streamflow for stations
  waterGap_streamflow = data.table::fread(watergap_raw_path)
  gauge_names = colnames(waterGap_streamflow)
  # defining the empty lists 
  list_min_p3m = list()
  list_mean_p3m = list()
  list_min_p12m = list()
  list_mean_p12m = list()
  sd_list = list()
  cv_list = list()
  
  # looping over stations to compute the high resolution predictors
  for (i in seq_along(gauge_names)) {
    WaterGap_past_3_month_min = waterGap_streamflow %>% 
      add_date(.) %>%
      dplyr::select(date, Q = gauge_names[i]) %>%
      min_past_3_month(.)
    WaterGap_past_3_month_mean = waterGap_streamflow %>% 
      add_date(.) %>%
      dplyr::select(date, Q = gauge_names[i]) %>%
      mean_past_3_month(.)
    WaterGap_past_12_month_min = waterGap_streamflow %>% 
      add_date(.) %>%
      dplyr::select(date, Q = gauge_names[i]) %>%
      min_past_12_month(.)
    WaterGap_past_12_month_mean = waterGap_streamflow %>% 
      add_date(.) %>%
      dplyr::select(date, Q = gauge_names[i]) %>%
      mean_past_12_month(.)
    df_mean = waterGap_streamflow %>% 
      add_date(.) %>%
      dplyr::select(date, Q = gauge_names[i]) %>% 
      mean_mon(.)
    df_sd = waterGap_streamflow %>% 
      add_date(.) %>%
      dplyr::select(date, Q = gauge_names[i]) %>% 
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
    list_min_p3m[[gauge_names[i]]] = WaterGap_past_3_month_min
    list_mean_p3m[[gauge_names[i]]] = WaterGap_past_3_month_mean
    list_min_p12m[[gauge_names[i]]] = WaterGap_past_12_month_min
    list_mean_p12m[[gauge_names[i]]] = WaterGap_past_12_month_mean
    sd_list[[gauge_names[i]]] = sd_ts
    cv_list[[gauge_names[i]]] = cv_ts
    
    if (i %% 100 == 0) {
      percentage_complete = round((i / length(gauge_names)) * 100, digits = 2)
      cat(paste("Iteration", i, "completed.", "Progress:", percentage_complete, "%"))
    }
  }
  # convert the lists to df
  cat("the lists of HR predictors are been converting to datatables.\n")
  
  # Streamflow of the month 
  df_debi = lapply(seq_along(list_min_p3m), function(i) {
    list_min_p3m[[i]] %>% 
      dplyr::filter(date >= start_date) %>% 
      dplyr::select(one_of("Q"))
  }) %>% 
    do.call("cbind",.) %>% 
    `colnames<-`(gauge_names) %>% 
    setDT() %>% 
    .[,dates :=seq.Date(start_date, end_date, "month")] %>% 
    melt(id.vars='dates', variable.name='gaugeid',
         value.name = "Q") %>% 
    .[,gaugeid :=as.character(gaugeid)]
  # minimum of the past 3 months
  df_min_p3m = lapply(seq_along(list_min_p3m), function(i) {
    list_min_p3m[[i]] %>% 
      dplyr::filter(date >= start_date) %>% 
      dplyr::select(one_of("min_p3m"))
  }) %>% 
    do.call("cbind",.) %>% 
    `colnames<-`(gauge_names) %>% 
    setDT() %>% 
    .[,dates :=seq.Date(start_date, end_date, "month")] %>% 
    melt(id.vars='dates', variable.name='gaugeid',
         value.name = "min_p3m") %>% 
    .[,gaugeid :=as.character(gaugeid)]
  # minimum of the past 12 months
  df_min_p12m = lapply(seq_along(list_min_p12m), function(i) {
    list_min_p12m[[i]] %>% 
      dplyr::filter(date >= start_date) %>% 
      dplyr::select(one_of("min_p12m"))
  }) %>% 
    do.call("cbind",.) %>% 
    `colnames<-`(gauge_names) %>% 
    setDT() %>% 
    .[,dates :=seq.Date(start_date, end_date, "month")] %>% 
    melt(id.vars='dates', variable.name='gaugeid',
         value.name = "min_p12m") %>% 
    .[,gaugeid :=as.character(gaugeid)]
  
  # mean of the past 3 months
  df_mean_p3m = lapply(seq_along(list_mean_p3m), function(i) {
    list_mean_p3m[[i]] %>% 
      dplyr::filter(date >= start_date) %>% 
      dplyr::select(one_of("mean_p3m"))
  }) %>% 
    do.call("cbind",.) %>% 
    `colnames<-`(gauge_names) %>% 
    setDT() %>% 
    .[,dates :=seq.Date(start_date, end_date, "month")] %>% 
    melt(id.vars='dates', variable.name='gaugeid',
         value.name = "mean_p3m") %>% 
    .[,gaugeid :=as.character(gaugeid)]
  # mean of the past 12 months
  df_mean_p12m = lapply(seq_along(list_mean_p12m), function(i) {
    list_mean_p12m[[i]] %>% 
      dplyr::filter(date >= start_date) %>% 
      dplyr::select(one_of("mean_p12m"))
  }) %>% 
    do.call("cbind",.) %>% 
    `colnames<-`(gauge_names) %>% 
    setDT() %>% 
    .[,dates :=seq.Date(start_date, end_date, "month")] %>% 
    melt(id.vars='dates', variable.name='gaugeid',
         value.name = "mean_p12m") %>% 
    .[,gaugeid :=as.character(gaugeid)]
  
  cat("the convertion of standard deviation and coffecient",
      "of varation have been left.")
  # sd
  df_sd = lapply(seq_along(sd_list), function(i) {
    sd_list[[i]] %>% 
      dplyr::filter(date >= start_date) %>% 
      dplyr::select(one_of("sd"))
  }) %>% 
    do.call("cbind",.) %>% 
    `colnames<-`(gauge_names) %>% 
    setDT() %>% 
    .[,dates :=seq.Date(start_date, end_date, "month")] %>% 
    melt(id.vars='dates', variable.name='gaugeid',
         value.name = "sd") %>% 
    .[,gaugeid :=as.character(gaugeid)]
  # cv
  df_cv = lapply(seq_along(cv_list), function(i) {
    cv_list[[i]] %>% 
      dplyr::filter(date >= start_date) %>% 
      dplyr::select(one_of("cv"))
  }) %>% 
    do.call("cbind",.) %>% 
    `colnames<-`(gauge_names) %>% 
    setDT() %>% 
    .[,dates :=seq.Date(start_date, end_date, "month")] %>% 
    melt(id.vars='dates', variable.name='gaugeid',
         value.name = "cv") %>% 
    .[,gaugeid :=as.character(gaugeid)]
  
  # combine all the predictors into a final list as output of the func.
  out <- cbind(df_debi, df_min_p3m[,.(min_p3m)], df_min_p12m[,.(min_p12m)],
               df_mean_p3m[, .(mean_p3m)], df_mean_p12m[, .(mean_p12m)],
               df_sd[, .(sd)], df_cv[,.(cv)])
  
  return(out)
}
### Low Resolution predictors -----
extract_LR_predictors = function(path, gauges_shp,
                                 variable_name = NULL,
                                 out_dir = NULL){
  
  filename <- list.files(path = path,
                         pattern = "\\.tif$",
                         full.names = TRUE) %>% 
    gtools::mixedsort()
  
  WaterGap_mat = matrix(NA, nrow = length(filename),
                        ncol = dim(gauges_shp)[1])
  
  for (i in seq_along(filename)) {
    
    r = raster(filename[i])
    rasValue = raster::extract(r, gauges_shp)
    WaterGap_mat[i,] = rasValue
    cat("The extraction of the LR vatiable for the gauging stations in",
        basename(filename[i]), "file was done.\n")
  }
  
  colnames(WaterGap_mat) = gauges_shp$gauge_d
  
  if (!file.exists(out_dir)){
    dir.create(out_dir, recursive = TRUE)
  }
  
  WaterGap_mat %>%
    data.table::as.data.table() %>%
    data.table::fwrite(., file = file.path(out_dir,
                                           paste0(variable_name, "_LR.csv")))
  
  
  return(WaterGap_mat)
}
### static predictors predictors -----
compute_static_predictors <- function(in_dt_path, gauge_shp){
  
  start_date <<- start_date
  end_date <<- end_date
  
  # Load the table of reaches and gauge shapefile
  data_net = data.table::fread(in_dt_path)
  
  # Pull the DRYvER ids and their corresponding gauge names
  shp_dryver_id = gauge_shp %>% pull(DRYVER_RIV)
  gauge_name = gauge_shp %>% pull(gauge_d)
  
  # Define an internal function to extend each predictor
  
  extend_pred <- function(data_net, shp_dryver_id, var_name,
                          gauge_name = gauge_name, number_iter = 468){
    
    start_date <<- start_date
    end_date <<- end_date
    
    if (var_name != "ai"){
      lapply(seq_along(shp_dryver_id), function(i){
        data_net[DRYVER_RIVID == shp_dryver_id[i]][[var_name]] %>%  rep(., number_iter)
      }) %>% 
        do.call("cbind",.) %>% 
        as.data.table() %>% 
        `colnames<-`(gauge_name) %>% 
        .[,dates :=seq.Date(start_date, end_date, "month")] %>%
        melt(id.vars='dates', variable.name='gaugeid',
             value.name = var_name) %>% 
        .[,gaugeid :=as.character(gaugeid)]
    } else {
      
      lapply(seq_along(shp_dryver_id), function(i){
        data_net[,.(DRYVER_RIVID, ai_ix_m01_uav, ai_ix_m02_uav, ai_ix_m03_uav, ai_ix_m04_uav, ai_ix_m05_uav,
                    ai_ix_m06_uav, ai_ix_m07_uav, ai_ix_m08_uav, ai_ix_m09_uav, ai_ix_m10_uav,
                    ai_ix_m11_uav, ai_ix_m12_uav)] %>% 
          .[DRYVER_RIVID == shp_dryver_id[i]] %>%  .[,2:13] %>% unlist() %>% as.vector() %>%  rep(., number_iter)
      }) %>% 
        do.call("cbind",.) %>% 
        as.data.table() %>% 
        `colnames<-`(gauge_name) %>% 
        .[,dates :=seq.Date(start_date, end_date, "month")] %>%
        melt(id.vars='dates', variable.name='gaugeid',
             value.name = var_name) %>% 
        .[,gaugeid :=as.character(gaugeid)]
    }
  }
  
  # Apply the internal function to all the static predictors and
  # convert them to the long format data.table
  ai_extended <- extend_pred(data_net = data_net, shp_dryver_id = shp_dryver_id,
                             var_name = 'ai', gauge_name = gauge_name, number_iter = 39)
  
  slope_extended <- extend_pred(data_net, shp_dryver_id = shp_dryver_id,
                                var_name = 'slope', gauge_name = gauge_name)
  
  glacier_fraction_extended <- extend_pred(data_net,shp_dryver_id = shp_dryver_id,
                                           var_name = 'glacier_fraction', gauge_name = gauge_name)
  
  drainage_area_extended <-  extend_pred(data_net, shp_dryver_id = shp_dryver_id,
                                         var_name = 'drainage_area', gauge_name = gauge_name)
  land_cover_extended <- extend_pred(data_net, shp_dryver_id = shp_dryver_id,
                                     var_name = 'land_cover', gauge_name = gauge_name) 
  
  pot_nat_vegetation_extended <- extend_pred(data_net, shp_dryver_id = shp_dryver_id,
                                             var_name = 'pot_nat_vegetation', gauge_name = gauge_name) 
  
  karst_fraction_extended <- extend_pred(data_net, shp_dryver_id = shp_dryver_id,
                                         var_name = 'karst_fraction', gauge_name = gauge_name) 
  karst_status_extended <- extend_pred(data_net, shp_dryver_id = shp_dryver_id,
                                       var_name = 'karst_status', gauge_name = gauge_name) 
  lka_pc_use_extended <- extend_pred(data_net, shp_dryver_id = shp_dryver_id,
                                     var_name = 'lka_pc_use', gauge_name = gauge_name) 
  ppd_pk_cav_extended <- extend_pred(data_net, shp_dryver_id = shp_dryver_id,
                                     var_name = 'ppd_pk_cav', gauge_name = gauge_name) 
  ppd_pk_uav_extended <- extend_pred(data_net, shp_dryver_id = shp_dryver_id,
                                     var_name = 'ppd_pk_uav', gauge_name = gauge_name) 
  ire_pc_cse_extended <- extend_pred(data_net, shp_dryver_id = shp_dryver_id,
                                     var_name = 'ire_pc_cse', gauge_name = gauge_name) 
  ire_pc_use_extended <- extend_pred(data_net, shp_dryver_id = shp_dryver_id,
                                     var_name = 'ire_pc_use', gauge_name = gauge_name) 
  dor_pc_pva_extended <- extend_pred(data_net, shp_dryver_id = shp_dryver_id,
                                     var_name = 'dor_pc_pva', gauge_name = gauge_name) 
  
  # Extract the coordinate of the stations to analyz spatial-temporal cv
  gauge_coordinates <- gauge_shp %>% sf::st_coordinates() %>% 
    as.data.table()
  
  x_coordinate <- gauge_coordinates[,X] %>% rep(.,468) %>%
    matrix(., ncol = 468) %>%
    t() %>% as.data.table() %>% setnames(gauge_name) %>% 
    .[,dates :=seq.Date(start_date, end_date, "month")] %>%
    melt(id.vars='dates', variable.name='gaugeid',
         value.name = "X")
  
  y_coordinate <- gauge_coordinates[,Y] %>% rep(.,468) %>%
    matrix(., ncol = 468) %>%
    t() %>% as.data.table() %>% setnames(gauge_name) %>% 
    .[,dates :=seq.Date(start_date, end_date, "month")] %>%
    melt(id.vars='dates', variable.name='gaugeid',
         value.name = "Y")
  
  function_output <- cbind(ai_extended, slope_extended[,.(slope)], glacier_fraction_extended[,.(glacier_fraction)],
                           drainage_area_extended[,.(drainage_area)], land_cover_extended[,.(land_cover)],
                           pot_nat_vegetation_extended[,.(pot_nat_vegetation)], karst_fraction_extended[,.(karst_fraction)],
                           karst_status_extended[,.(karst_status)], lka_pc_use_extended[,.(lka_pc_use)],
                           ppd_pk_cav_extended[,.(ppd_pk_cav)], ppd_pk_uav_extended[,.(ppd_pk_uav)], ire_pc_cse_extended[,.(ire_pc_cse)],
                           ire_pc_use_extended[,.(ire_pc_use)], dor_pc_pva_extended[,.(dor_pc_pva)],
                           x_coordinate[,.(X)], y_coordinate[,.(Y)])
  
  return(function_output)
}

combined_predictors <- function(in_target, in_hr_pred, path_lr_pred, in_static_pred){
  
  start_date <<- start_date
  end_date <<- end_date
  # Compute the number of no-flow days per each month for stations
  
  in_target[, year_month := format(dates, "%Y-%m")]
  in_target[value < 0.001, value := 0]
  target_dt <- in_target[, .(target = sum(value == 0)), by = .(gaugeid, year_month)][
    , dates := as.Date(paste0(year_month, "-01"))][
      ,.(gaugeid, dates, target)
      ]
  keep_gaugesid <- target_dt[,gaugeid] %>% unique()
  # Import and compute the LR predictors 
  lr_pre_files <- list.files(path = path_lr_pred,
                             pattern = ".csv",
                             full.names = TRUE)
  out <- lapply(seq_along(lr_pre_files), function(i){
    data.table::fread(lr_pre_files[i]) %>% 
      .[,dates :=seq.Date(start_date, end_date, "month")] %>% 
      melt(id.vars='dates', variable.name='gaugeid',
           value.name = paste0("V", i)) %>% 
      .[,gaugeid :=as.character(gaugeid)]
  })
  
  lr_pred_dt <- cbind(out[[1]], out[[2]][,.(V2)]) %>% 
    setnames(c("V1", "V2"),
             c("gwr_to_runoff_ratio", "wet_days")) %>% 
    .[gaugeid %in% keep_gaugesid]
  
  # Combine the predictors and target into one table
  combined_pred <- cbind(in_hr_pred[gaugeid %in% keep_gaugesid], lr_pred_dt[,.(gwr_to_runoff_ratio, wet_days)],
                         in_static_pred[,-c(1,2)])
  
  # Divid the predictors with unites m3/sec to mm/sec
  cols_to_divid <- c("Q", "min_p3m", "min_p12m", "mean_p3m", "mean_p12m")
  combined_pred[,(cols_to_divid) := lapply(.SD, function(x) x/drainage_area),
                .SDcols = cols_to_divid]
  
  function_output <- data.table::merge.data.table(target_dt, combined_pred,
                                                  by = c("gaugeid", "dates"), all = TRUE) %>% 
    data.table::setorder(., "gaugeid") %>% 
    .[complete.cases(.)]
  
  return(function_output)
  
}


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
    po_classweights <- mlr3pipelines::po("classweights", minor_weight = imbalance_ratio)
    #Create graph learners so that oversampling happens systematically upstream of all training
    lrns[['lrn_ranger_oversampled']] <- mlr3pipelines::GraphLearner$new(
      po_oversampled %>>% lrns[['lrn_ranger']])
    lrns[['lrn_ranger_classweights']] <- mlr3pipelines::GraphLearner$new(
      po_classweights %>>% lrns[['lrn_ranger']])
  }
  
  return(lrns)
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
                   lower = floor(0.1*nfeatures),
                   upper = floor(0.5*nfeatures)), #Half number of features
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
                                          pvalue = TRUE, pvalue_permutn = 25) {
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
                                     pvalue = TRUE, pvalue_permutn = 25) {
  
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


# ggvimp ------------
ggvimp <- function(in_rftuned, in_predvars, varnum = 17, spatial_rsp=FALSE) {
  rsmp_res <- get_outerrsmp(in_rftuned, spatial_rsp=spatial_rsp)
  
  #Get variable importance and format them
  varimp_basic <- weighted_vimportance_nestedrf(rfresamp = rsmp_res,
                                                pvalue = FALSE) %>%
    merge(., in_predvars, by.x='varnames', by.y = "varname") %>%
    .[, `:=`(varnames = factor(varnames, varnames[order(-imp_wmean)]),
             Category = factor(Category,
                               levels = c('Climate', 'Hydrology', 'Landcover', 'Man-made',
                                          'Physiography', 'Soils & Geology'))
    )] %>%
    setorder(-imp_wmean)
  # varimp_basic$varnames <- c("drainage area", "slope", "Q_sd_12_s", "Q_s", "Q_min_p12_s", "P_to_PET_ratio",
  #                            "Natural vegetation", "runoff_dvar_s", "Q_mean_p3_s","gwr_to_runoff_ratio",
  #                            "Q_cv_12_s","Q_min_p3_s","Karst_extent_s","Q_mean_p12_s",
  #                            "Land cover",  "wet days", "karst status") %>% as.factor(.)
  outp <- varimp_basic[1:varnum] %>% 
    mutate(varnames = fct_reorder(varnames, desc(imp_wmean))) %>% 
    ggplot(. ,aes(x=varnames,
                  color =Category, fill=Category)) +
    geom_bar(aes(y=imp_wmean), stat = 'identity', alpha=0.7) +
    geom_errorbar(aes(ymin=imp_wmean-imp_wsd, ymax=imp_wmean+imp_wsd)) +
    scale_x_discrete(labels = function(x) {
      stringr::str_wrap(tolower(x), width = 27)
    },
    limits=rev) +
    scale_fill_manual(values=c('#fdb462','#80b1d3','#b3de69','#bc80bd','#696868', '#696868'),
                      drop=FALSE) +
    scale_color_manual(values=c('#fdb462','#80b1d3','#b3de69','#bc80bd','#696868', '#696868'),
                       drop=FALSE) +
    theme_classic() +
    theme(axis.text.x = element_text(size=8),
          # axis.title.x  = element_blank(),
          legend.text = element_text(size=12),
          legend.title = element_text(size=14),
          legend.position = c(0.80, 0.5)
    ) +
    scale_y_continuous(expand=c(0,0), position = 'left') +
    labs(y = "Variable importance (actual impurity reduction)")+
    coord_flip(ylim=c(0, max(varimp_basic[, max(imp_wmean+imp_wsd)+10], 100))) 
  #Plot 'em
  
  
  return(outp)
}

predname_df <- function(in_task, feat_name_vec = NULL, category = NULL){
  
  if (inherits(in_task, "Task")) {
    feat_name <- in_task$feature_names
  } else
    feat_name <- feat_name_vec
  
  if (is.null(category)) {
    cat <- c("Hydrology", "Climate", "Hydrology", "Man-made",  "Physiography","Hydrology",
             "Man-made", "Man-made", "Soils & Geology", "Soils & Geology",
             "Landcover",  "Landcover", "Hydrology", "Hydrology", "Hydrology", "Hydrology",
             "Landcover", "Man-made", "Man-made", "Hydrology", "Physiography", "Climate")
  } else
    cat <- category
  
  predvars <- cbind(feat_name, cat) %>%
    `colnames<-`(c("varname", "Category"))
  
  return(predvars)
  
}
# extract_pd_nestedrf -----------
extract_pd_nestedrf <- function(learner_id=1, in_rftuned, datdf,
                                selcols, nvariate, grid_resolution) {
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
  
  # ngridvec <- c(ngrid, ngrid)
  
  #Make dataset of all combinations of selected column names, two at a time
  if (nvariate == 1) {
    
    pdcomb <- lapply(selcols, function(i) {
      print(i)
      
      pdout <- pdp::partial(in_fit,
                            pred.var = c(i), 
                            train = datdf,
                            grid.resolution = grid_resolution) %>% 
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
      
      pdout <- pdp::partial(in_fit,
                            pred.var = c(i, j), 
                            train = datdf,
                            grid.resolution = grid_resolution) %>%
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
ggpartialdep <- function (in_rftuned, in_predvars, colnums, grid_resolution=10, nodupli=F,
                          nvariate = 1, parallel=T, spatial_rsp=FALSE) {
  
  #Get outer resampling of interest
  rsmp_res <- get_outerrsmp(in_rftuned, spatial_rsp=spatial_rsp)
  
  #Get partial dependence across all folds and repeats
  nlearners <-with(rsmp_res$resampling$param_set$values, folds*repeats)
  datdf <- as.data.frame(rsmp_res$task$data()) #This may be shortened
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
                                      grid_resolution = grid_resolution,
                                      future.scheduling = structure(TRUE,ordering = "random"),
                                      future.packages = c("data.table","pdp","ranger"))
    
  } else {
    print(paste("Computing partial dependence iteratively across", nlearners,
                "CV folds"))
    pd <- lapply(seq_len(nlearners),
                 extract_pd_nestedrf,
                 in_rftuned = rsmp_res,
                 datdf = datdf,
                 selcols = selcols,
                 nvariate = nvariate,
                 grid_resolution = grid_resolution)
  }
  
  #Get weighted mean
  varvec <- paste0('var', 1:nvariate)
  valvec <- paste0('value', 1:nvariate)
  
  pdformat <- do.call(rbind, pd) %>%
    setDT %>%
    .[, list(mean1 = weighted.mean(`yhat`, classif.bacc)),
      by= c(varvec, valvec)] %>%
    .[, variables := var1] %>%
    merge(., in_predvars, by.x='var1', by.y='varname')
  
  datdf2 <- as.data.table(datdf)[, target := as.numeric(as.character(target))]
  
  if (nvariate ==1) {
    tileplots_l <- pdformat[,list(list(ggplotGrob(
      ggplot(.SD, aes(x=value1, y=mean1)) +
        geom_line() +
        geom_rug(data=datdf2,
                 aes_string(x=eval(var1),y='target'),
                 alpha=1/3) +
        scale_y_continuous(name='Partial dependence (probability of intermittency)',
                           limits= c(min(mean1)-0.01, max(mean1)+0.01),  #c(0.25, 0.425),
                           expand=c(0,0))+
        scale_x_continuous(name=stringr::str_wrap(eval(variables), width = 30)) +
        theme_classic() +
        theme(text = element_text(size=12),
              axis.title.y = element_blank())
    ))), by=.(var1)]
    
  } else if (nvariate == 2) {
    pdformat[, variables := paste(var1, var2)]
    
    vargrid <- t(combn(1:length(selcols), 2))
    #leglims <- pdformat[, c(min(mean1), max(mean1))]
    
    #Iterate over every pair of variables
    
    tileplots_l <- pdformat[,list(list(ggplotGrob(
      ggplot(.SD, aes(x=value1, y=value2)) +
        geom_tile(aes(fill = mean1)) +
        scale_fill_distiller(palette='Grey') +
        geom_jitter(data=datdf,
                    aes_string(color='intermittent_o1800', x=eval(var1),y=eval(var2)),
                    alpha=1/3) +
        scale_color_manual(values=c('#0F9FD6','#ff9b52')) +
        labs(x=stringr::str_wrap(in_predvars[varcode==eval(var1), varname],
                                 width = 20),
             y=stringr::str_wrap(in_predvars[varcode==eval(var2), varname],
                                 width = 20)) +
        theme_bw() +
        theme(text = element_text(size=12))
    )))
    , by=.(var1, var2)]
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
      left = 'Partial dependence (probability of intermittency)')))
  })
  return(tileplots_multipl)
  }

