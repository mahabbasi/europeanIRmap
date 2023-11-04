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
    c('sm_1180', 'smires', 'record is too incomplete', 'months with zeros', '2011-01-01', NA),
    c('sm_1098', 'smires', 'no-data for the reach- static predictors', 'remove entire record', NA, NA)
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
  
  # drainage_area_extended <-  extend_pred(data_net, shp_dryver_id = shp_dryver_id,
  #                                        var_name = 'drainage_area', gauge_name = gauge_name)
  
  drainage_area_extended <-  lapply(seq_along(shp_dryver_id), function(i){
    gauge_name_vec <- gauge_name[i]
    gauge_shp %>%
      sf::st_drop_geometry() %>%
      as.data.table() %>%
      .[DRYVER_RIV == shp_dryver_id[i]] %>% 
      .[gauge_d == gauge_name_vec, up_mdfd] %>%  rep(., 468)
  }) %>% 
    do.call("cbind",.) %>% 
    as.data.table() %>% 
    `colnames<-`(gauge_name) %>% 
    .[,dates :=seq.Date(start_date, end_date, "month")] %>%
    melt(id.vars='dates', variable.name='gaugeid',
         value.name = 'drainage_area') %>% 
    .[,gaugeid :=as.character(gaugeid)]
  
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
