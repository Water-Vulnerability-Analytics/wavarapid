#' Calculate Runoff Exposure Indicator (Standardized Version)
#' 
#' Calculates runoff vulnerability using LOCA2 climate model data.
#' First checks if processed files exist, if not runs NetCDF processing.
#' 
#' @param water_supplies sf object with water supply locations
#' @param runoff_data character vector of processed loca2 mwbm files
#' @param save logical, whether to save output
#' @param output_dir character, directory path for saving files
#' @param verbose logical, whether to show progress messages
#' 
#' @return tibble with water supply ID, park info, and runoff exposure indicators
#' 
calc_runoff_exposure <- function(water_supplies,
                                 runoff_data = "data/rapid/inputs/loca2_mwbm_runoff/",
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
    message("Starting runoff exposure calculation...")
    start_time <- Sys.time()
  }
  
  # pull files
  files <- list.files(runoff_data, full.names = TRUE)
  
  # Check if we have any files to process
  if (length(files) == 0) {
    stop("No processed files found. Check your input data and file paths.")
  }
  
  if (verbose) {
    message("Processing ", length(files), " runoff files")
    
    # Create progress bar for file processing
    pb <- progress::progress_bar$new(
      format = "Processing runoff files [:bar] :percent (:current/:total) ETA: :eta",
      total = length(files),
      clear = FALSE,
      width = 80
    )
  }
  
  # Process the extracted data files with progress bar
  mean_data <- map(files, ~{
    if (verbose) pb$tick()
    process_raw_grid_data(.x, water_supplies = water_supplies)
  })
  
  final_data <- mean_data %>% bind_rows()
  
  # Separate historical and future data
  historical <- final_data %>% 
    filter(period == "training") %>%
    rename(runoff_historical = mean_runoff)
  
  future <- final_data %>% 
    filter(period != "training") %>%
    rename(runoff_future = mean_runoff) 
  
  # Calculate runoff change
  runoff_change <- left_join(future, historical, by = c("wsd_source_id", "park_unit")) %>%
    select(park_unit, wsd_source_id, model = model.x, ssp = period.x, run = run.x, 
           runoff_historical, period_historical = time.x, runoff_future, period_future = time.y) %>%
    mutate(runoff_change = runoff_future - runoff_historical,
           percent_change = (runoff_change / runoff_historical) * 100)
  
  # Calculate percentile indicators using standardized naming
  p10_values <- runoff_change %>%
    group_by(park_unit, wsd_source_id) %>%
    mutate(p10 = quantile(percent_change, 0.10, na.rm = TRUE),
           diff_from_p10 = abs(percent_change - p10)) %>%
    filter(diff_from_p10 == min(diff_from_p10, na.rm = TRUE)) %>%
    slice(1) %>%  # Take first if multiple ties
    ungroup() %>%
    # Keep original column names
    transmute(
      wsd_source_id,
      park_unit,
      raw_runoff_change_p10 = percent_change
    )
  
  p50_values <- runoff_change %>%
    group_by(park_unit, wsd_source_id) %>%
    mutate(p50 = median(percent_change, na.rm = TRUE),
           diff_from_p50 = abs(percent_change - p50)) %>%
    filter(diff_from_p50 == min(diff_from_p50, na.rm = TRUE)) %>%
    slice(1) %>%  # Take first if multiple ties
    ungroup() %>%
    transmute(
      wsd_source_id,
      park_unit,
      raw_runoff_change_p50 = percent_change
    )
  
  p90_values <- runoff_change %>%
    group_by(park_unit, wsd_source_id) %>%
    mutate(p90 = quantile(percent_change, 0.90, na.rm = TRUE),
           diff_from_p90 = abs(percent_change - p90)) %>%
    filter(diff_from_p90 == min(diff_from_p90, na.rm = TRUE)) %>%
    slice(1) %>%  # Take first if multiple ties
    ungroup() %>%
    transmute(
      wsd_source_id,
      park_unit,
      raw_runoff_change_p90 = percent_change
    )
  
  # Calculate fraction of models showing decrease
  model_fractions <- runoff_change %>%
    group_by(park_unit, wsd_source_id) %>%
    summarize(
      raw_models_showing_decrease_pct = 100 * sum(percent_change < 0, na.rm = TRUE) / n(),
      .groups = "drop"
    )
  
  # Combine all indicators using standardized joins
  final_results <- p10_values %>%
    left_join(p50_values, by = c("wsd_source_id", "park_unit")) %>%
    left_join(p90_values, by = c("wsd_source_id", "park_unit")) %>%
    left_join(model_fractions, by = c("wsd_source_id", "park_unit")) %>%
    # Add park names from water_supplies
    left_join(
      water_supplies %>% 
        st_drop_geometry() %>% 
        select(wsd_source_id, park_unit, park_name),
      by = c("wsd_source_id", "park_unit")
    ) %>%
    # Add standardized indicators using helper function
    add_standardized_indicators("raw_runoff_change") %>%
    add_standardized_indicators("raw_models_showing_decrease") %>%
    # Reorder columns following standardized pattern
    select(wsd_source_id, park_unit, park_name, everything())
  
  # Save results
  if (save) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    filename <- file.path(output_dir, paste0("runoff_exposure_", Sys.Date(), ".csv"))
    write_csv(final_results, filename)
    
    if (verbose) {
      message("Results saved to: ", filename)
    }
  }
  
  # Print timing
  if (verbose) {
    end_time <- Sys.time()
    total_time <- end_time - start_time
    message("\nRunoff exposure processing completed in ", 
            round(total_time, 2), " ", attr(total_time, "units"))
  }
  
  return(final_results)
}