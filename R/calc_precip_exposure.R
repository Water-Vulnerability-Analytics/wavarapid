#' @title Precipitation exposure: projected change in annual precipitation
#'
#' @description
#' Calculates the change in mean annual precipitation at each water supply
#' between a historical and a future period, from LOCA2 monthly
#' precipitation. Monthly precipitation is averaged over each site's grid
#' cells and summed to annual totals. Only files whose years (from the file
#' name) overlap the two periods are read. For each model run (model, scenario and
#' ensemble member) the change is the future mean minus the historical mean,
#' as a percentage of the historical mean.
#'
#' The ensemble is then summarized per site: the runs closest to the 10th,
#' 50th and 90th percentiles of the percent change, and the percentage of runs
#' that show a decrease.
#'
#' @param water_supplies sf object. The water supplies, with a
#'   `wsd_source_id` column, e.g. from [wavacore::get_wsd()]. Only their IDs
#'   are used, to pick rows from the processed files, so they must be the IDs
#'   the files were made with. Required.
#' @param precip_data String. The folder of processed precipitation parquet
#'   files from [process_loca2_mwbm()]. Stops with an error if no files are
#'   left after the `hindcast` and year filters. Default:
#'   `rapid/inputs/loca2_precipitation_water_supplies` in the shared data
#'   folder,
#'   see [wavacore::wava_data_path()].
#' @param hindcast Logical. Where the historical period comes from.
#'   * `TRUE`: each model's own historical run, so the change is measured
#'     within each model. Historical years after the historical run ends
#'     (2014) come from the same model's scenario run. Training and Livneh
#'     (`livneh_unsplit`) files are not read.
#'   * `FALSE`: observed (Livneh) precipitation, the same baseline for every
#'     model, from the `livneh_unsplit` files in `precip_data`; the future
#'     comes from the scenario files. Model values before 2019 are dropped.
#'     Needs Livneh files in `precip_data`.
#'
#'   Default: `TRUE`.
#' @param historical Character or numeric vector of 2 years, e.g.
#'   `c(1980, 2010)`. The first and last year of the historical period.
#'   Required.
#' @param future Character or numeric vector of 2 years, e.g.
#'   `c(2035, 2065)`. The first and last year of the future period. Required.
#' @param water_year Logical. How months are grouped into years.
#'   * `TRUE`: water years, October to September, labelled by the year they
#'     end in (water year 2010 is October 2009 to September 2010).
#'   * `FALSE`: calendar years.
#'
#'   Default: `TRUE`.
#' @param save Logical. Whether to write the results to a CSV file.
#'   * `TRUE`: writes `precipitation_exposure_<date>.csv`, named with today's date,
#'     to `output_dir`, creating the folder if needed. Stops with an error if
#'     `output_dir` is `NULL`.
#'   * `FALSE`: writes nothing.
#'
#'   Default: `FALSE`.
#' @param output_dir String or `NULL`. The folder for the CSV file.
#'   * A folder path: where the file is written when `save = TRUE`. Ignored
#'     when `save = FALSE`.
#'   * `NULL`: only allowed with `save = FALSE`.
#'
#'   Default: `NULL`.
#' @param verbose Logical. Whether to print progress.
#'   * `TRUE`: shows a progress bar over the files and messages with the run
#'     time and the path of the saved file.
#'   * `FALSE`: runs silently.
#'
#'   Default: `TRUE`.
#'
#' @return A tibble, one row per water supply, with `wsd_source_id`,
#'   `park_unit` and:
#'   * `raw_precipitation_change_p10`, `_p50`, `_p90`: the percent change in
#'     mean annual precipitation of the model run closest to that percentile;
#'   * `raw_precip_models_showing_decrease`: the percentage of model runs
#'     with a decrease;
#'   * `norm_*` and `rank_*` versions of both from
#'     [add_standardized_indicators()]. The precipitation change is inverted, so the
#'     largest decrease gets the highest score.
#' @export
#' @seealso [calc_precip_sensitivity()] for the historical trend;
#'   [calc_precip_timing_exposure()] for the change in timing;
#'   [process_loca2_mwbm()], which makes the input files;
#'   [calc_vulnerability_index()].
#' @examplesIf interactive() && dir.exists(wavacore::wava_data_path("rapid", "inputs", "loca2_precipitation_water_supplies"))
#' # Reads the processed precipitation files for both periods (several GB),
#' # so it takes several minutes even for one park
#' ws <- wavacore::get_wsd(park = "ZION")
#' precip_exp <- calc_precip_exposure(ws,
#'                                    historical = c(1980, 2010),
#'                                    future = c(2035, 2065))
#' precip_exp

calc_precip_exposure <- function(water_supplies,
                                 precip_data = wavacore::wava_data_path("rapid", "inputs", "loca2_precipitation_water_supplies/"),
                                 hindcast = TRUE,
                                 historical,
                                 future,
                                 water_year = TRUE,
                                 save = FALSE,
                                 output_dir = NULL,
                                 verbose = TRUE) {
  
  # Input validation
  stopifnot(
    "water_supplies must be an sf object" = inherits(water_supplies, "sf"),
    "save is TRUE but output_dir not provided" = !save || !is.null(output_dir)
  )
  
  # Start timing
  if (verbose) {
    message("Starting precipitation exposure calculation...")
    start_time <- Sys.time()
  }
  
  # Extract year boundaries from date arguments
  hist_start_year <- as.integer(historical[1])
  hist_end_year <- as.integer(historical[2])
  future_start_year <- as.integer(future[1])
  future_end_year <- as.integer(future[2])
  
  if(hindcast == FALSE){
    # if we want to use training data... select associated files
    files <- list.files(precip_data, full.names = TRUE) %>%
      grep("2015|2045|2075|livneh_unsplit", ., value = TRUE)
    
  } else if(hindcast == TRUE){
    # if we want to use hindcast climate data... 
    files <- list.files(precip_data, full.names = TRUE) %>%
      grep("training|livneh_unsplit", .,  invert = TRUE, value = TRUE)
  } 
  
  files_dt <- data.table::data.table(file = files)
  
  # extract years
  files_dt[, `:=`(
    year_single = as.integer(stringr::str_extract(file, "(?<!\\d)(\\d{4})(?=\\.nc)")),
    year_start = as.integer(stringr::str_extract(file, "(\\d{4})(?=-\\d{4})")),
    year_end = as.integer(stringr::str_extract(file, "(?<=-)(\\d{4})(?=\\.)"))
  )]
  
  # fill single-year files into start/end
  files_dt[!is.na(year_single), `:=`(
    year_start = year_single,
    year_end = year_single
  )]
  
  if(water_year == FALSE){
    # keep files overlapping either period
    files_subset <- files_dt[
      (year_end >= hist_start_year & year_start <= hist_end_year) |
        (year_end >= future_start_year  & year_start <= future_end_year)
    ]
  } else {
    # keep files overlapping either period (with wyear accounted for)
    files_subset <- files_dt[
      (year_end >= hist_start_year - 1 & year_start <= hist_end_year) |
        (year_end >= future_start_year - 1  & year_start <= future_end_year)
    ]
  }
  
  files <- files_subset$file
  
  # Check if we have any files to process
  if (length(files) == 0) {
    stop("No processed files found. Check your input data and file paths.")
  }
  
  if (verbose) {
    message("Processing ", length(files), " precipitation files")
    
    # Create progress bar for file processing
    pb <- progress::progress_bar$new(
      format = "Processing precipitation files [:bar] :percent (:current/:total) ETA: :eta",
      total = length(files),
      clear = FALSE,
      width = 80)
  }
  
  # function that employs data.table to get annual averages from the raw tabular
  # data for every water supply, calculates the mean precipitation across neighborhood
  # grids and times
  process_raw_grid_data <- function(f) {
    
    raw <- arrow::read_parquet(f) %>%
      dplyr::filter(wsd_source_id %in% water_supplies$wsd_source_id)
    
    annual_points <- data.table::as.data.table(raw)[
      ,
      year := if (water_year) {
        lubridate::year(date) + data.table::fifelse(lubridate::month(date) >= 10, 1L, 0L)
      } else {
        lubridate::year(date)
      }
    ]
    
    # when using observed data, get rid of the modeled data 
    # that overlaps with observed record
    if (!grepl("livneh_unsplit", f) && hindcast == FALSE) {
      annual_points <- annual_points[year >= 2019]
    }
    
    annual_points <- annual_points[
      # Use the extracted year boundaries from function arguments
      (year >= hist_start_year & year <= hist_end_year) |
        (year >= future_start_year & year <= future_end_year)
    ][
      # Average whichever neighborhood cells are available at each time step.
      , .(site_precipitation = if (any(is.finite(value))) {
        mean(value[is.finite(value)])
      } else {
        NA_real_
      }),
      by = .(wsd_source_id, date, year, park_unit, nc_file)
    ][
      # Sum the resulting site series; an entirely missing year remains NA.
      , .(total_precipitation = if (any(is.finite(site_precipitation))) {
        sum(site_precipitation[is.finite(site_precipitation)])
      } else {
        NA_real_
      }),
      by = .(year, wsd_source_id, park_unit, nc_file)
    ]
    
    annual_points
  }
  
  # Process the extracted data files with progress bar
  mean_data <- purrr::map(files, ~{
    if (verbose) pb$tick(tokens = list(what = "Processing files"))
    process_raw_grid_data(.x)
  })
  
  combined_data <- mean_data %>% dplyr::bind_rows()
  
  data.table::setDT(combined_data)  
  
  suppressMessages({suppressWarnings({
    # fix specific filename pattern from Jay's data, if using
    combined_data[
      nc_file == "Jay_Alder_pr_MWBM_LOCA2_training_195001-201412.nc",
      nc_file := "Jay_Alder_pr_MWBM_LOCA2.historical.training.1950-2100.nc"
    ]
    
    # parse filename components
    combined_data[
      , c("var", "model", "ssp", "run", "nc", "misc1", "misc2", "misc3" ) := data.table::tstrsplit(nc_file, "\\.")
    ]
    
    # ssp populate with "historical" for hindcast = FALSE (as in, no ssp)
    combined_data[
      !ssp %in% c("ssp245", "ssp370", "ssp585", "historical"),
      ssp := dplyr::case_when(
        grepl("livneh_unsplit", nc_file, ignore.case = TRUE) ~ "historical",
        TRUE ~ stringr::str_extract(nc_file, "ssp\\d+")
      )
    ]
    
    # Split run/time if present
    combined_data[
      !is.na(run), c("run", "time") := data.table::tstrsplit(run, "_", fixed = TRUE)
    ]
  })})
  
  # Drop unnecessary columns
  combined_data[, c("var", "nc") := NULL]
  
  if(hindcast == TRUE){
    
    # distinct (model, scenario) pairs that exist in the past and future
    ssps <- combined_data %>%
      dplyr::filter(!ssp %in% c("1950-2014", "historical")) %>%
      dplyr::distinct(model, ssp, run)
    
    precipitation_change_sub <- vector("list", length = nrow(ssps))
    
    for(i in 1:nrow(ssps)){
      
      ssp_sub <- ssps[i,]
      
      hindcast_historical <- combined_data %>%
        dplyr::filter(model == ssp_sub$model,
               ssp == "historical",
               run == ssp_sub$run) %>%
        dplyr::filter(year >= hist_start_year & year <= hist_end_year)
      
      # if our historical record is actually considered "future" data...
      hindcast_future <- combined_data %>%
        dplyr::filter(model == ssp_sub$model,
               ssp == ssp_sub$ssp,
               run == ssp_sub$run) %>%
        dplyr::filter(year >= hist_start_year & year <= hist_end_year,
               !year %in% unique(hindcast_historical$year))
      
      historical_complete <- hindcast_historical %>%
        dplyr::bind_rows(hindcast_future) %>%
        dplyr::group_by(wsd_source_id, park_unit) %>%
        dplyr::summarize(precipitation_historical = mean(total_precipitation, na.rm = TRUE),
                  range_historical = paste0(min(year), "-", max(year)),
                  .groups = "drop")
      
      future_complete <- combined_data %>% 
        dplyr::filter(model == ssp_sub$model,
               ssp == ssp_sub$ssp,
               run == ssp_sub$run) %>%
        dplyr::filter(year >= future_start_year & year <= future_end_year) %>%
        dplyr::group_by(wsd_source_id, park_unit) %>%
        dplyr::summarize(precipitation_future = mean(total_precipitation, na.rm = TRUE),
                  period_future = paste0(min(year), "-", max(year)),
                  .groups = "drop")
      
      precipitation_change_sub[[i]] <- future_complete %>%
        dplyr::inner_join(., historical_complete, by = c("wsd_source_id", "park_unit")) %>%
        dplyr::mutate(model = ssp_sub$model,
               ssp = ssp_sub$ssp,
               run = ssp_sub$run,
               precipitation_change = precipitation_future - precipitation_historical,
               percent_change = (precipitation_change / precipitation_historical) * 100)
    }
    
    precipitation_change <- dplyr::bind_rows(precipitation_change_sub)
    
  } else if(hindcast == FALSE){
    
    # distinct (model, scenario) pairs that exist in the past and future
    ssps <- combined_data %>%
      dplyr::filter(ssp != "training" & ssp != "historical") %>%
      dplyr::distinct(model, ssp, run)
    
    training_historical <- combined_data %>%
      dplyr::filter(ssp == "training" | ssp == "historical") %>%
      dplyr::group_by(wsd_source_id, park_unit) %>%
      dplyr::filter(year >= hist_start_year & year <= hist_end_year) %>%
      dplyr::summarize(precipitation_historical = mean(total_precipitation, na.rm = TRUE),
                range_historical = paste0(min(year), "-", max(year)),
                .groups = "drop")
    
    precipitation_change_sub <- vector("list", length = nrow(ssps))
    
    for(i in 1:nrow(ssps)){
      suppressMessages({suppressWarnings({
        ssp_sub <- ssps[i,]
        
        future_complete <- combined_data %>% 
          dplyr::filter(model == ssp_sub$model,
                 ssp == ssp_sub$ssp,
                 run == ssp_sub$run) %>%
          dplyr::filter(year >= future_start_year & year <= future_end_year) %>%
          dplyr::group_by(wsd_source_id, park_unit) %>%
          dplyr::summarize(precipitation_future = mean(total_precipitation, na.rm = TRUE),
                    period_future = paste0(min(year), "-", max(year)),
                    .groups = "drop")
        
        precipitation_change_sub[[i]] <- future_complete %>%
          dplyr::inner_join(., training_historical, by = c("wsd_source_id", "park_unit")) %>%
          dplyr::mutate(model = ssp_sub$model,
                 ssp = ssp_sub$ssp,
                 run = ssp_sub$run,
                 precipitation_change = precipitation_future - precipitation_historical,
                 percent_change = (precipitation_change / precipitation_historical) * 100)
        
      })})
      
    }
    
    precipitation_change <- dplyr::bind_rows(precipitation_change_sub)
    
  }
  
  # calculate percentile indicators:
  suppressMessages({suppressWarnings({
    
    p10_values <- precipitation_change %>%
      dplyr::group_by(park_unit, wsd_source_id) %>%
      dplyr::mutate(p10 = stats::quantile(percent_change, 0.10, na.rm = TRUE),
             diff_from_p10 = abs(percent_change - p10)) %>%
      dplyr::filter(diff_from_p10 == min(diff_from_p10, na.rm = TRUE)) %>%
      dplyr::slice(1) %>%  # take first if ties
      dplyr::ungroup() %>%
      dplyr::select(wsd_source_id,
             park_unit,
             raw_precipitation_change_p10 = percent_change)
    
    p50_values <- precipitation_change %>%
      dplyr::group_by(park_unit, wsd_source_id) %>%
      dplyr::mutate(p50 = stats::median(percent_change, na.rm = TRUE),
             diff_from_p50 = abs(percent_change - p50)) %>%
      dplyr::filter(diff_from_p50 == min(diff_from_p50, na.rm = TRUE)) %>%
      dplyr::slice(1) %>%  # take first if ties
      dplyr::ungroup() %>%
      dplyr::select(wsd_source_id,
             park_unit,
             raw_precipitation_change_p50 = percent_change)
    
    p90_values <- precipitation_change %>%
      dplyr::group_by(park_unit, wsd_source_id) %>%
      dplyr::mutate(p90 = stats::quantile(percent_change, 0.90, na.rm = TRUE),
             diff_from_p90 = abs(percent_change - p90)) %>%
      dplyr::filter(diff_from_p90 == min(diff_from_p90, na.rm = TRUE)) %>%
      dplyr::slice(1) %>%  # take first if ties
      dplyr::ungroup() %>%
      dplyr::select(wsd_source_id,
             park_unit,
             raw_precipitation_change_p90 = percent_change)
    
    # calculate fraction of models showing decrease
    model_fractions <- precipitation_change %>%
      dplyr::group_by(park_unit, wsd_source_id) %>%
      dplyr::summarize(raw_precip_models_showing_decrease = 100 * sum(percent_change < 0, na.rm = TRUE) / dplyr::n(),
                .groups = "drop")
    
  })})
  
  # Combine all indicators using standardized joins
  final_results <- p10_values %>%
    dplyr::left_join(p50_values, by = c("wsd_source_id", "park_unit")) %>%
    dplyr::left_join(p90_values, by = c("wsd_source_id", "park_unit")) %>%
    dplyr::left_join(model_fractions, by = c("wsd_source_id", "park_unit")) %>%
    # add standardized indicators using helper function
    add_standardized_indicators("raw_precipitation_change", invert = TRUE) %>%
    add_standardized_indicators("raw_precip_models_showing_decrease") %>%
    # reorder columns following standardized pattern
    dplyr::select(wsd_source_id, park_unit, dplyr::everything())
  
  # Save results
  if (save) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    filename <- file.path(output_dir, paste0("precipitation_exposure_", Sys.Date(), ".csv"))
    readr::write_csv(final_results, filename)
    
    if (verbose) {
      message("Results saved to: ", filename)
    }
  }
  
  # Print timing
  if (verbose) {
    end_time <- Sys.time()
    total_time <- end_time - start_time
    message("\nprecipitation exposure processing completed in ",
            round(total_time, 2), " ", attr(total_time, "units"))
  }
  
  return(final_results)
}
