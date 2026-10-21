#' Calculate Snow Exposure Indicator
#' Calculates peak SWE change vulnerability using LOCA2 climate model data and the MWBM.
#'
#' @param water_supplies sf object with water supply locations
#' @param swe_data character vector of processed loca2 mwbm files
#' @param hindcast logical, if TRUE uses hindcast climate data, if FALSE uses "real" climate data 
#' @param historical Character vector of length 2 giving the start and end years
#'   (in `"YYYY"` format) that defines the "historical" period.  
#' @param future Character vector of length 2 giving the start and end years
#'   (in `"YYYY"` format) that defines the "future" period.  
#' @param water_year T/F, if TRUE, uses water-year bounds (Oct 1 - Sep 30) before annual aggregation.
#' @param save logical, whether to save output
#' @param output_dir character, directory path for saving files
#' @param verbose logical, whether to show progress messages
#'
#' @return tibble with water supply ID, park info, and SWE exposure indicators
#'
#' @export
calc_swe_exposure <- function(water_supplies,
                              swe_data = wavacore::wava_data_path("rapid", "inputs", "loca2_swe_water_supplies/"),
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
    message("Starting SWE exposure calculation...")
    start_time <- Sys.time()
  }
  
  # Extract year boundaries from date arguments
  hist_start_year <- as.integer(historical[1])
  hist_end_year <- as.integer(historical[2])
  future_start_year <- as.integer(future[1])
  future_end_year <- as.integer(future[2])
  
  if(hindcast == FALSE){
    # if we want to use training data... select associated files
    files <- list.files(swe_data, full.names = TRUE) %>%
      grep("2015-2100|training", ., value = TRUE)
    
  } else if(hindcast == TRUE){
    # if we want to use hindcast climate data... 
    files <- list.files(swe_data, full.names = TRUE) %>%
      grep("training", .,  invert = TRUE, value = TRUE)
  } 
  
  # Check if we have any files to process
  if (length(files) == 0) {
    stop("No processed files found. Check your input data and file paths.")
  }
  
  if (verbose) {
    message("Processing ", length(files), " SWE files")
    
    # Create progress bar for file processing
    pb <- progress::progress_bar$new(
      format = "Processing SWE files [:bar] :percent (:current/:total) ETA: :eta",
      total = length(files),
      clear = FALSE,
      width = 80)
  }
  
  # function that employs data.table to get annual avg from the raw tabular
  # data for every water supply, calculates the mean SWE across neighborhood
  # grids and times
  process_raw_grid_data <- function(f) {
    
    raw <- arrow::read_parquet(f) %>%
      dplyr::filter(wsd_source_id %in% water_supplies$wsd_source_id)
    
    # grabs year (or water year if water_year = TRUE)
    annual_points <- data.table::as.data.table(raw)[
      ,      year := if (water_year) {
        lubridate::year(date) + data.table::fifelse(lubridate::month(date) >= 10, 1L, 0L)
      } else {
        lubridate::year(date)
      }
    ][
      # Use the extracted year boundaries from function arguments
      (year >= hist_start_year & year <= hist_end_year) | 
        (year >= future_start_year & year <= future_end_year)
    ][
      # mean SWE across grid cells per month
      , .(mean_swe = if (any(is.finite(value))) {
        mean(value[is.finite(value)])
      } else {
        NA_real_
      }),
      by = .(wsd_source_id, date, year, park_unit, nc_file)
    ][
      # average annual max per water source/model
      , .(max_swe = if (any(is.finite(mean_swe))) {
        max(mean_swe[is.finite(mean_swe)])
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
  
  suppressMessages({suppressWarnings({
    # fix specific filename pattern from Jay's data
    combined_data[
      nc_file == "Jay_Alder_snow_MWBM_LOCA2_training_195001-201412.nc",
      nc_file := "Jay_Alder_snow_MWBM_LOCA2.historical.training.1950-2100.nc"
    ]
    
    # parse filename components
    combined_data[
      , c("var", "model", "ssp", "run", "nc") := data.table::tstrsplit(nc_file, "\\.")
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
      dplyr::filter(time != "1950-2014") %>%
      dplyr::distinct(model, ssp, run)
    
    swe_change_sub <- vector("list", length = nrow(ssps))
    
    for(i in 1:nrow(ssps)){
      
      ssp_sub <- ssps[i,]
      
      hindcast_historical <- combined_data %>%
        dplyr::filter(model == ssp_sub$model,
               ssp == "historical",
               run == ssp_sub$run) %>%
        dplyr::filter(year >= hist_start_year & year <= hist_end_year)
      
      hindcast_future <- combined_data %>%
        dplyr::filter(model == ssp_sub$model,
               ssp == ssp_sub$ssp,
               run == ssp_sub$run) %>%
        dplyr::filter(year >= hist_start_year & year <= hist_end_year)
      
      historical_complete <- hindcast_historical %>%
        dplyr::bind_rows(hindcast_future) %>%
        dplyr::group_by(wsd_source_id, park_unit) %>%
        dplyr::summarize(swe_historical = mean(max_swe, na.rm = TRUE),
                  range_historical = paste0(min(year), "-", max(year)),
                  .groups = "drop")
      
      future_complete <- combined_data %>% 
        dplyr::filter(model == ssp_sub$model,
               ssp == ssp_sub$ssp,
               run == ssp_sub$run) %>%
        dplyr::filter(year >= future_start_year & year <= future_end_year) %>%
        dplyr::group_by(wsd_source_id, park_unit) %>%
        dplyr::summarize(swe_future = mean(max_swe, na.rm = TRUE),
                  period_future = paste0(min(year), "-", max(year)),
                  .groups = "drop")
      
      swe_change_sub[[i]] <- future_complete %>%
        dplyr::inner_join(., historical_complete, by = c("wsd_source_id", "park_unit")) %>%
        dplyr::mutate(model = ssp_sub$model,
               ssp = ssp_sub$ssp,
               run = ssp_sub$run,
               swe_change = swe_future - swe_historical,
               percent_change = (swe_change / swe_historical) * 100)
    }
    
    swe_change <- dplyr::bind_rows(swe_change_sub)
    
  } else if(hindcast == FALSE){
    
    # distinct (model, scenario) pairs that exist in the past and future
    ssps <- combined_data %>%
      dplyr::filter(ssp != "training") %>%
      dplyr::distinct(model, ssp, run)
    
    swe_change_sub <- vector("list", length = nrow(ssps))
    
    for(i in 1:nrow(ssps)){
      suppressMessages({suppressWarnings({
        ssp_sub <- ssps[i,]
        
        training_historical <- combined_data %>%
          dplyr::filter(ssp == "training") %>%
          dplyr::filter(year >= hist_start_year & year <= hist_end_year)
        
        model_historical <- combined_data %>%
          dplyr::filter(model == ssp_sub$model,
                 ssp == ssp_sub$ssp,
                 run == ssp_sub$run) %>%
          dplyr::filter(year >= hist_start_year & year <= hist_end_year)
        
        historical_complete <- training_historical %>%
          dplyr::bind_rows(model_historical) %>%
          dplyr::group_by(wsd_source_id, park_unit) %>%
          dplyr::summarize(swe_historical = mean(max_swe, na.rm = TRUE),
                    range_historical = paste0(min(year), "-", max(year)),
                    .groups = "drop")
        
        future_complete <- combined_data %>% 
          dplyr::filter(model == ssp_sub$model,
                 ssp == ssp_sub$ssp,
                 run == ssp_sub$run) %>%
          dplyr::filter(year >= future_start_year & year <= future_end_year) %>%
          dplyr::group_by(wsd_source_id, park_unit) %>%
          dplyr::summarize(swe_future = mean(max_swe, na.rm = TRUE),
                    period_future = paste0(min(year), "-", max(year)),
                    .groups = "drop")
        
        swe_change_sub[[i]] <- future_complete %>%
          dplyr::inner_join(., historical_complete, by = c("wsd_source_id", "park_unit")) %>%
          dplyr::mutate(model = ssp_sub$model,
                 ssp = ssp_sub$ssp,
                 run = ssp_sub$run,
                 swe_change = swe_future - swe_historical,
                 percent_change = (swe_change / swe_historical) * 100)
        
      })})
      
    }
    
    swe_change <- dplyr::bind_rows(swe_change_sub)
    
  }
  
  # calculate percentile indicators:
  suppressMessages({suppressWarnings({
    
    p10_values <- swe_change %>%
      dplyr::group_by(park_unit, wsd_source_id) %>%
      dplyr::mutate(p10 = stats::quantile(percent_change, 0.10, na.rm = TRUE),
             diff_from_p10 = abs(percent_change - p10)) %>%
      dplyr::filter(diff_from_p10 == min(diff_from_p10, na.rm = TRUE)) %>%
      dplyr::slice(1) %>%  # take first if ties
      dplyr::ungroup() %>%
      dplyr::select(wsd_source_id,
             park_unit,
             raw_swe_change_p10 = percent_change)
    
    p50_values <- swe_change %>%
      dplyr::group_by(park_unit, wsd_source_id) %>%
      dplyr::mutate(p50 = stats::median(percent_change, na.rm = TRUE),
             diff_from_p50 = abs(percent_change - p50)) %>%
      dplyr::filter(diff_from_p50 == min(diff_from_p50, na.rm = TRUE)) %>%
      dplyr::slice(1) %>%  # take first if ties
      dplyr::ungroup() %>%
      dplyr::select(wsd_source_id,
             park_unit,
             raw_swe_change_p50 = percent_change)
    
    p90_values <- swe_change %>%
      dplyr::group_by(park_unit, wsd_source_id) %>%
      dplyr::mutate(p90 = stats::quantile(percent_change, 0.90, na.rm = TRUE),
             diff_from_p90 = abs(percent_change - p90)) %>%
      dplyr::filter(diff_from_p90 == min(diff_from_p90, na.rm = TRUE)) %>%
      dplyr::slice(1) %>%  # take first if ties
      dplyr::ungroup() %>%
      dplyr::select(wsd_source_id,
             park_unit,
             raw_swe_change_p90 = percent_change)
    
    # calculate fraction of models showing decrease
    model_fractions <- swe_change %>%
      dplyr::group_by(park_unit, wsd_source_id) %>%
      dplyr::summarize(raw_swe_models_showing_decrease = 100 * sum(percent_change < 0, na.rm = TRUE) / dplyr::n(),
                .groups = "drop")
    
  })})
  
  # Combine all indicators using standardized joins
  final_results <- p10_values %>%
    dplyr::left_join(p50_values, by = c("wsd_source_id", "park_unit")) %>%
    dplyr::left_join(p90_values, by = c("wsd_source_id", "park_unit")) %>%
    dplyr::left_join(model_fractions, by = c("wsd_source_id", "park_unit")) %>%
    # add standardized indicators using helper function
    add_standardized_indicators("raw_swe_change", invert = TRUE) %>%
    add_standardized_indicators("raw_swe_models_showing_decrease") %>%
    # reorder columns following standardized pattern
    dplyr::select(wsd_source_id, park_unit, dplyr::everything())
  
  # Save results
  if (save) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    filename <- file.path(output_dir, paste0("swe_exposure_", Sys.Date(), ".csv"))
    readr::write_csv(final_results, filename)
    
    if (verbose) {
      message("Results saved to: ", filename)
    }
  }
  
  # Print timing
  if (verbose) {
    end_time <- Sys.time()
    total_time <- end_time - start_time
    message("\nSWE exposure processing completed in ",
            round(total_time, 2), " ", attr(total_time, "units"))
  }
  
  return(final_results)
}
