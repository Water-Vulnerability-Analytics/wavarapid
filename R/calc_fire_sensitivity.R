#' Calculate Fire Sensitivity Indicator
#' 
#' Calculates fire sensitivity using Wildfire Hazard Potential (WHP) data
#' from Dillon et al. 2023. Computes area-weighted mean WHP within water
#' supply areas.
#' 
#' @param water_supplies sf object with water supply locations
#' @param whp_raster terra raster of wildfire hazard potential
#' @param save logical, whether to save output
#' @param output_dir character, directory path for saving files
#' @param verbose logical, whether to show progress messages
#' 
#' @return tibble with water supply ID, park info, and fire sensitivity indicators
#' 
calc_fire_sensitivity <- function(water_supplies,
                                  whp_raster,
                                  save = FALSE,
                                  output_dir = NULL,
                                  verbose = TRUE) {
  
  # Input validation
  stopifnot(
    "water_supplies must be an sf object" = inherits(water_supplies, "sf"),
    "whp_raster must be a terra raster" = inherits(whp_raster, "SpatRaster"),
    "save is TRUE but output_dir not provided" = !save || !is.null(output_dir)
  )
  
  # Start timing
  start_time <- Sys.time()
  
  # Create progress bar
  if (verbose) {
    pb <- progress::progress_bar$new(
      format = "Processing fire sensitivity [:bar] :percent (:current/:total) ETA: :eta",
      total = nrow(water_supplies),
      clear = FALSE,
      width = 80
    )
  }
  
  # Process each water supply
  results <- map_dfr(seq_len(nrow(water_supplies)), function(i) {
    if (verbose) pb$tick()
    
    # Calculate area-weighted mean WHP for water supply buffer
    ws_whp <- terra::zonal(
      whp_raster,
      terra::vect(st_transform(water_supplies[i, ], crs(whp_raster))),
      fun = "mean",
      weights = TRUE,
      na.rm = TRUE
    )
    
    # Create results tibble
    tibble(
      wsd_source_id = water_supplies[i, ]$wsd_source_id,
      park_unit = water_supplies[i, ]$park_unit,
      park_name = water_supplies[i, ]$park_name,
      raw_wildfire_hazard_mean = pull(ws_whp)
    )
  })
  
  # Calculate standardized indicators using helper function  
  final_results <- results %>%
    add_standardized_indicators("raw_wildfire_hazard")
  
  # Save results
  if (save) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    filename <- file.path(output_dir, paste0("fire_sensitivity_", Sys.Date(), ".csv"))
    write_csv(final_results, filename)
    
    if (verbose) {
      message("Results saved to: ", filename)
    }
  }
  
  # Print timing
  if (verbose) {
    end_time <- Sys.time()
    total_time <- end_time - start_time
    message("\nFire sensitivity processing completed in ", 
            round(total_time, 2), " ", attr(total_time, "units"))
  }
  
  return(final_results)
}