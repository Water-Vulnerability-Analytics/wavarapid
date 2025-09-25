#' Calculate Fire Exposure Indicator
#' 
#' Calculates future fire probability changes for water supplies using
#' data from Gao et al. 2021. Computes area-weighted future fire probability
#' for water supply pyromes across multiple climate scenarios.
#' 
#' @param water_supplies sf object with water supply locations
#' @param fire_prob_data fire probability data from Gao et al. 2021
#' @param pyrome_data sf object of CONUS pyromes
#' @param save logical, whether to save output
#' @param output_dir character, directory path for saving files
#' @param verbose logical, whether to show progress messages
#' 
#' @return tibble with water supply ID, park info, and fire exposure indicators
#' 
calc_fire_exposure <- function(water_supplies,
                               fire_prob_data,
                               pyrome_data,
                               save = FALSE,
                               output_dir = NULL,
                               verbose = TRUE) {
  
  # Input validation
  stopifnot(
    "water_supplies must be an sf object" = inherits(water_supplies, "sf"),
    "save is TRUE but output_dir not provided" = !save || !is.null(output_dir)
  )
  
  # Start timing
  start_time <- Sys.time()
  
  # convert both objects to terra for faster filtering
  pyrome_data <- vect(pyrome_data)
  water_supplies <- vect(water_supplies)
  
  # Create progress bar
  if (verbose) {
    pb <- progress::progress_bar$new(
      format = "Processing fire exposure [:bar] :percent (:current/:total) ETA: :eta",
      total = nrow(water_supplies),
      clear = FALSE,
      width = 80
    )
  }
  
  # Process each water supply
  results <- map_dfr(seq_len(nrow(water_supplies)), function(i) {
    if (verbose) pb$tick()
    
    # Get water supply pyrome
    # ws_pyrome <- pyrome_data %>% 
    #   st_filter(water_supplies[i, ])
    ws_pyrome <- pyrome_data[relate(pyrome_data, water_supplies[i,], "intersects"), ]
    
    # Extract fire probability data for this pyrome
    fp_data <- fire_prob_data %>% 
      filter(PYROME == ws_pyrome$PYROME) %>% 
      select(delta_fp_p10, delta_fp_p50, delta_fp_p90) %>% 
      mutate(
        wsd_source_id = water_supplies[i, ]$wsd_source_id,
        park_unit = water_supplies[i, ]$park_unit,
        park_name = water_supplies[i, ]$park_name
      ) %>% 
      select(wsd_source_id, park_unit, park_name, everything())
    
    return(fp_data)
  })
  
  # Calculate standardized indicators
  final_results <- results %>%
    # Raw values (renamed for consistency)
    rename(
      raw_fire_prob_change_p10 = delta_fp_p10,
      raw_fire_prob_change_p50 = delta_fp_p50,
      raw_fire_prob_change_p90 = delta_fp_p90
    ) %>%
    # Use helper function to add normalized and ranked columns
    add_standardized_indicators("raw_fire_prob_change")
  
  # Save results
  if (save) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    filename <- file.path(output_dir, paste0("fire_exposure_", Sys.Date(), ".csv"))
    write_csv(final_results, filename)
    
    if (verbose) {
      message("Results saved to: ", filename)
    }
  }
  
  # Print timing
  if (verbose) {
    end_time <- Sys.time()
    total_time <- end_time - start_time
    message("\nFire exposure processing completed in ", 
            round(total_time, 2), " ", attr(total_time, "units"))
  }
  
  return(final_results)
}