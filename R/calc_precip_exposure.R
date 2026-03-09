#' Calculate Precipitation Exposure Indicator (Standardized Version)
#'
#' Calculates precipitation vulnerability using LOCA2 climate model data and the MWBM.
#'
#' @param water_supplies sf object with water supply locations
#' @param precip_data character vector of processed loca2 mwbm files
#' @param hindcast logical, if TRUE uses hindcast climate data, if FALSE uses "real" climate data 
#' @param historical Character vector of length 2 giving the start and end dates
#'   (in `"YYYY-MM-DD"` format) that defines the "historical" period.  
#' @param future Character vector of length 2 giving the start and end dates
#'   (in `"YYYY-MM-DD"` format) that defines the "future" period.  
#' @param save logical, whether to save output
#' @param output_dir character, directory path for saving files
#' @param verbose logical, whether to show progress messages
#'
#' @return tibble with water supply ID, park info, and precipitation exposure indicators
#'
calc_precip_exposure <- function(water_supplies,
                                 precip_data = "data/rapid/inputs/loca2_precipitation_water_supplies/",
                                 hindcast = TRUE,
                                 historical = c("1981-01-01", "2010-12-31"),
                                 future = c("2035-01-01", "2065-12-31"),
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
  hist_start_year <- as.integer(format(as.Date(historical[1]), "%Y"))
  hist_end_year <- as.integer(format(as.Date(historical[2]), "%Y"))
  future_start_year <- as.integer(format(as.Date(future[1]), "%Y"))
  future_end_year <- as.integer(format(as.Date(future[2]), "%Y"))
  
  if(hindcast == FALSE){
    # if we want to use training data... select associated files
    files <- list.files(precip_data, full.names = TRUE) %>%
      grep("2015-2100|training", ., value = TRUE)
    
  } else if(hindcast == TRUE){
    # if we want to use hindcast climate data... 
    files <- list.files(precip_data, full.names = TRUE) %>%
      grep("training", .,  invert = TRUE, value = TRUE)
  } 
  
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
    
    raw <- read_parquet(f) %>%
      dplyr::filter(wsd_source_id %in% water_supplies$wsd_source_id)
    
    annual_points <- as.data.table(raw)[
      , year := as.integer(format(date, "%Y"))
    ][
      # Use the extracted year boundaries from function arguments
      (year >= hist_start_year & year <= hist_end_year) | 
        (year >= future_start_year & year <= future_end_year)
    ][
      # annual totals per grid cell/model
      , .(grid_precipitation = sum(value, na.rm = TRUE)),
      by = .(wsd_source_id, year, park_unit, nc_file, cell_lon, cell_lat)
    ][
      # average annual totals per water source/model
      , .(total_precipitation = mean(grid_precipitation, na.rm = TRUE)),
      by = .(year, wsd_source_id, park_unit, nc_file)
    ]
    
  }
  
  # Process the extracted data files with progress bar
  mean_data <- map(files, ~{
    if (verbose) pb$tick(tokens = list(what = "Processing files"))
    process_raw_grid_data(.x)
  })
  
  combined_data <- mean_data %>% bind_rows()
  
  suppressMessages({suppressWarnings({
    # fix specific filename pattern from Jay's data
    combined_data[
      nc_file == "Jay_Alder_pr_MWBM_LOCA2_training_195001-201412.nc",
      nc_file := "Jay_Alder_pr_MWBM_LOCA2.historical.training.1950-2100.nc"
    ]
    
    # parse filename components
    combined_data[
      , c("var", "model", "ssp", "run", "nc", "misc1", "misc2", "misc3" ) := tstrsplit(nc_file, "\\.")
    ]
    
    # Split run/time if present
    combined_data[
      !is.na(run), c("run", "time") := tstrsplit(run, "_", fixed = TRUE)
    ]
  })})
  
  # Drop unnecessary columns
  combined_data[, c("var", "nc") := NULL]
  
  if(hindcast == TRUE){
    
    # distinct (model, scenario) pairs that exist in the past and future
    ssps <- combined_data %>%
      filter(time != "1950-2014") %>%
      distinct(model, ssp, run)
    
    precipitation_change_sub <- vector("list", length = nrow(ssps))
    
    for(i in 1:nrow(ssps)){
      
      ssp_sub <- ssps[i,]
      
      hindcast_historical <- combined_data %>%
        filter(model == ssp_sub$model,
               ssp == "historical",
               run == ssp_sub$run) %>%
        filter(year >= hist_start_year & year <= hist_end_year)
      
      hindcast_future <- combined_data %>%
        filter(model == ssp_sub$model,
               ssp == ssp_sub$ssp,
               run == ssp_sub$run) %>%
        filter(year >= hist_start_year & year <= hist_end_year)
      
      historical <- hindcast_historical %>%
        bind_rows(hindcast_future) %>%
        group_by(wsd_source_id, park_unit) %>%
        summarize(precipitation_historical = mean(total_precipitation, na.rm = TRUE),
                  range_historical = paste0(min(year), "-", max(year)),
                  .groups = "drop")
      
      future <- combined_data %>% 
        filter(model == ssp_sub$model,
               ssp == ssp_sub$ssp,
               run == ssp_sub$run) %>%
        filter(year >= future_start_year & year <= future_end_year) %>%
        group_by(wsd_source_id, park_unit) %>%
        summarize(precipitation_future = mean(total_precipitation, na.rm = TRUE),
                  period_future = paste0(min(year), "-", max(year)),
                  .groups = "drop")
      
      precipitation_change_sub[[i]] <- future %>%
        inner_join(., historical, by = c("wsd_source_id", "park_unit")) %>%
        mutate(model = ssp_sub$model,
               ssp = ssp_sub$ssp,
               run = ssp_sub$run,
               precipitation_change = precipitation_future - precipitation_historical,
               percent_change = (precipitation_change / precipitation_historical) * 100)
    }
    
    precipitation_change <- bind_rows(precipitation_change_sub)
    
  } else if(hindcast == FALSE){
    
    # distinct (model, scenario) pairs that exist in the past and future
    ssps <- combined_data %>%
      filter(ssp != "training") %>%
      distinct(model, ssp, run)
    
    precipitation_change_sub <- vector("list", length = nrow(ssps))
    
    for(i in 1:nrow(ssps)){
      suppressMessages({suppressWarnings({
        ssp_sub <- ssps[i,]
        
        training_historical <- combined_data %>%
          filter(ssp == "training") %>%
          filter(year >= hist_start_year & year <= hist_end_year)
        
        model_historical <- combined_data %>%
          filter(model == ssp_sub$model,
                 ssp == ssp_sub$ssp,
                 run == ssp_sub$run) %>%
          filter(year >= hist_start_year & year <= hist_end_year)
        
        historical <- training_historical %>%
          bind_rows(model_historical) %>%
          group_by(wsd_source_id, park_unit) %>%
          summarize(precipitation_historical = mean(total_precipitation, na.rm = TRUE),
                    range_historical = paste0(min(year), "-", max(year)),
                    .groups = "drop")
        
        future <- combined_data %>% 
          filter(model == ssp_sub$model,
                 ssp == ssp_sub$ssp,
                 run == ssp_sub$run) %>%
          filter(year >= future_start_year & year <= future_end_year) %>%
          group_by(wsd_source_id, park_unit) %>%
          summarize(precipitation_future = mean(total_precipitation, na.rm = TRUE),
                    period_future = paste0(min(year), "-", max(year)),
                    .groups = "drop")
        
        precipitation_change_sub[[i]] <- future %>%
          inner_join(., historical, by = c("wsd_source_id", "park_unit")) %>%
          mutate(model = ssp_sub$model,
                 ssp = ssp_sub$ssp,
                 run = ssp_sub$run,
                 precipitation_change = precipitation_future - precipitation_historical,
                 percent_change = (precipitation_change / precipitation_historical) * 100)
        
      })})
      
    }
    
    precipitation_change <- bind_rows(precipitation_change_sub)
    
  }
  
  # calculate percentile indicators:
  suppressMessages({suppressWarnings({
    
    p10_values <- precipitation_change %>%
      group_by(park_unit, wsd_source_id) %>%
      mutate(p10 = quantile(percent_change, 0.10, na.rm = TRUE),
             diff_from_p10 = abs(percent_change - p10)) %>%
      filter(diff_from_p10 == min(diff_from_p10, na.rm = TRUE)) %>%
      slice(1) %>%  # take first if ties
      ungroup() %>%
      select(wsd_source_id,
             park_unit,
             raw_precipitation_change_p10 = percent_change)
    
    p50_values <- precipitation_change %>%
      group_by(park_unit, wsd_source_id) %>%
      mutate(p50 = median(percent_change, na.rm = TRUE),
             diff_from_p50 = abs(percent_change - p50)) %>%
      filter(diff_from_p50 == min(diff_from_p50, na.rm = TRUE)) %>%
      slice(1) %>%  # take first if ties
      ungroup() %>%
      select(wsd_source_id,
             park_unit,
             raw_precipitation_change_p50 = percent_change)
    
    p90_values <- precipitation_change %>%
      group_by(park_unit, wsd_source_id) %>%
      mutate(p90 = quantile(percent_change, 0.90, na.rm = TRUE),
             diff_from_p90 = abs(percent_change - p90)) %>%
      filter(diff_from_p90 == min(diff_from_p90, na.rm = TRUE)) %>%
      slice(1) %>%  # take first if ties
      ungroup() %>%
      select(wsd_source_id,
             park_unit,
             raw_precipitation_change_p90 = percent_change)
    
    # calculate fraction of models showing decrease
    model_fractions <- precipitation_change %>%
      group_by(park_unit, wsd_source_id) %>%
      summarize(raw_models_showing_decrease_pct = 100 * sum(percent_change < 0, na.rm = TRUE) / n(),
                .groups = "drop")
    
  })})
  
  # Combine all indicators using standardized joins
  final_results <- p10_values %>%
    left_join(p50_values, by = c("wsd_source_id", "park_unit")) %>%
    left_join(p90_values, by = c("wsd_source_id", "park_unit")) %>%
    left_join(model_fractions, by = c("wsd_source_id", "park_unit")) %>%
    # add standardized indicators using helper function
    add_standardized_indicators("raw_precipitation_change", invert = TRUE) %>%
    add_standardized_indicators("raw_precip_models_showing_decrease") %>%
    # reorder columns following standardized pattern
    select(wsd_source_id, park_unit, everything())
  
  # Save results
  if (save) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    filename <- file.path(output_dir, paste0("precipitation_exposure_", Sys.Date(), ".csv"))
    write_csv(final_results, filename)
    
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